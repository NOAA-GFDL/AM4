module river_mod
!-----------------------------------------------------------------------
!                   GNU General Public License
!
! This program is free software; you can redistribute it and/or modify it and
! are expected to follow the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2 of
! the License, or (at your option) any later version.
!
! MOM is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.
! or see:   http://www.gnu.org/licenses/gpl.html
!-----------------------------------------------------------------------
! <CONTACT EMAIL="Kirsten.Findell@@noaa.gov"> Kirsten Findell </CONTACT>
! <CONTACT EMAIL="Zhi.Liang@@noaa.gov"> Zhi Liang </CONTACT>
! <NAMELIST NAME="river_nml">
! <DATA NAME="layout" TYPE="integer, dimension(2)">
!  Processor domain layout for river model. If layout(1)*layout(2) is not equal
!  to mpp_npes, the river model layout will be assigned the layout of land model
!  passed through river_init.
!  </DATA>
! <DATA NAME="do_rivers" TYPE="logical">
!   set true to run river model ( default is true). If FALSE, rivers are
!   essentially turned off to save computing time
!  </DATA>
! <DATA NAME="dt_slow" TYPE="real">
!   slow time step for river model. dt_slow must be integer multiplier of dt_fast passed
!   from land model ( land model time step).
!  </DATA>
! <DATA NAME="diag_freq" TYPE="integer">
!   Number of slow time steps between sending out diagnositics data(default is 1). Please
!   note that diagnostic output frequency ( specified in diag_table ) must be divided by
!   diag_freq*dt_slow.
!  </DATA>
! </NAMELIST>

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file
#endif

  use mpp_mod,             only : CLOCK_SUBCOMPONENT, CLOCK_ROUTINE
  use mpp_mod,             only : mpp_error, FATAL, WARNING, NOTE, stdout, stdlog
  use mpp_mod,             only : mpp_pe, mpp_chksum, mpp_max
  use mpp_mod,             only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, MPP_CLOCK_DETAILED
  use mpp_domains_mod,     only : domain2d, mpp_get_compute_domain, mpp_get_global_domain
  use mpp_domains_mod,     only : mpp_get_data_domain, mpp_update_domains, mpp_get_ntile_count
  use mpp_domains_mod,     only : domainUG, mpp_get_UG_compute_domain, mpp_pass_ug_to_sg
  use mpp_domains_mod,     only : mpp_pass_sg_to_ug
  use fms_mod,             only : check_nml_error, string
  use fms_mod,             only : close_file, file_exist, field_size, read_data, write_data
  use fms_mod,             only : field_exist, CLOCK_FLAG_DEFAULT
  use fms_io_mod,          only : get_mosaic_tile_file, get_instance_filename
  use fms_io_mod,          only : restart_file_type, register_restart_field, restore_state, save_restart, free_restart_type
  use diag_manager_mod,    only : diag_axis_init, register_diag_field, register_static_field, send_data, diag_field_add_attribute
  use time_manager_mod,    only : time_type, increment_time, get_time
  use data_override_mod,   only : data_override
  use river_type_mod,      only : river_type, Leo_Mad_trios, NO_RIVER_FLAG
  use river_physics_mod,   only : river_physics_step, river_physics_init, river_impedes_lake, &
                                  river_impedes_large_lake
  use constants_mod,       only : PI, RADIAN, tfreeze, DENS_H2O, hlf
  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT
  use land_tile_mod,       only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
  use land_data_mod,       only : land_data_type, lnd_sg, log_version, lnd
  use lake_tile_mod,       only : num_l
  use field_manager_mod, only: fm_field_name_len, fm_string_len, &
     fm_type_name_len, fm_path_name_len, fm_dump_list, fm_get_length, &
     fm_get_current_list, fm_loop_over_list, fm_change_list
  use land_io_mod,         only : new_land_io
  use fm_util_mod, only : fm_util_get_real, fm_util_get_logical, fm_util_get_string
  use tracer_manager_mod, only : NO_TRACER
  use table_printer_mod

  implicit none
  private

!--- version information ---------------------------------------------
character(len=*), parameter :: module_name = 'river_mod'
#include "../shared/version_variable.inc"

!--- public interface ------------------------------------------------
  public :: river_init, river_end, river_type, update_river, river_stock_pe
  public :: save_river_restart
  public :: river_tracers_init
  public :: num_river_tracers
  public :: river_tracer_index

!--- namelist interface ----------------------------------------------
  logical            :: do_rivers       = .TRUE.  ! if FALSE, rivers are essentially turned off to save computing time
  real               :: dt_slow
  integer            :: diag_freq       = 1       ! Number of slow time steps between sending out diagnositics data.
  logical            :: debug_river     = .FALSE.
  real               :: Somin           = 0.00005 ! There are 7 points with So = -9.999 but basinid > 0....
  real               :: outflowmean_min = 1.      ! temporary fix, should not allow zero in input file
  logical            :: land_area_called_cellarea = .false.
  logical            :: all_big_outlet_ctn0 = .false.

  real, dimension(3) :: ave_DHG_exp = (/0.49,0.33,0.18/)  ! (/B, F, M for avg of many rivers, 15Nov05/)
  real, dimension(3) :: ave_AAS_exp = (/0.19,0.39,0.42/)  ! (/b, f, m for avg of many rivers, 15Nov05/)
  real, dimension(3) :: ave_DHG_coef = (/4.62,0.26,0.82/) ! (/A, C, K for avg of many rivers, 15Nov05/)
  real               :: sinuosity = 1.3
  real               :: channel_tau = 86400*365.25*10     ! channel geometry reflects average flow over O(10 y)
  logical :: lake_area_bug = .FALSE. ! if set to true, reverts to buggy (quebec)
      ! behavior, where by mistake cell area was used instead of land area to
      ! compute the area of lakes.
  logical :: stop_on_mask_mismatch = .TRUE. ! if set to false, then the data mismatches (mmismatch
      ! of land and river masks, and discharges in pouints where there is no ocean) are reported,
      ! but don't cause the abort of the program.
  ! ZMS
  logical :: tracers_from_runoff = .false. ! if true, use runoff_c(:,:,num_phys+1:num_species)
          ! rather than source concentration and flux files

  namelist /river_nml/ dt_slow, diag_freq, debug_river,                      &
                       Somin, outflowmean_min, ave_DHG_exp, ave_AAS_exp,     &
                       ave_DHG_coef, do_rivers, sinuosity, channel_tau,      &
                       land_area_called_cellarea, all_big_outlet_ctn0,       &
                       lake_area_bug, stop_on_mask_mismatch,                 &
                       tracers_from_runoff

  character(len=128) :: river_src_file   = 'INPUT/river_data.nc'
  character(len=128) :: river_Omean_file = 'INPUT/river_Omean.nc'
!---------------------------------------------------------------------
  logical :: module_is_initialized = .FALSE.
  integer :: isc, iec, jsc, jec                         ! compute domain decomposition
  integer :: isd, ied, jsd, jed                         ! data domain decomposition
  integer :: lsc, lec                                   ! unstructure domain decomposition
  integer :: nlon, nlat                                 ! size of computational river grid
  integer :: num_lake_lev
  integer :: id_outflowmean, id_lake_depth_sill
  integer :: id_dx, id_basin, id_So, id_depth, id_width, id_vel
  integer :: id_LWSr, id_FWSr, id_HSr, id_meltr
  integer :: id_travel, id_elev, id_tocell
  integer :: maxtravel
  real    :: missing = -1.e8

  real,    parameter :: CONST_OMEAN = 80000
  real,    parameter :: epsln = 1.e-6
  real,    parameter :: sec_in_day = 86400.

  real     :: discharge_tol=0.0, clw=0.0, csw=0.0  ! will get these values from land model
  integer  :: i_river_ice, i_river_heat, i_river_DOC
  logical, allocatable, dimension(:,:) :: missing_rivers
  real,  allocatable, dimension(:,:)   :: discharge2ocean_next   ! store discharge value
  real,  allocatable, dimension(:,:,:) :: discharge2ocean_next_c ! store discharge value
  ! IDs of diag fields normalized per land area
  integer, allocatable, dimension(:)   :: id_infloc,  id_storage, id_stordis, id_inflow, &
        id_run_stor, id_outflow, id_removal, id_dis, id_lake_outflow
  ! IDs of diag fields normalized per cell area
  integer, allocatable, dimension(:)   :: id_infloc_c,  id_storage_c, id_stordis_c, id_inflow_c, &
        id_run_stor_c, id_outflow_c, id_removal_c, id_dis_c, id_lake_outflow_c
  integer :: id_dis_liq,  id_dis_ice,  id_dis_heat, id_dis_sink, id_dis_DOC, id_no_riv
  integer :: num_fast_calls
  integer :: slow_step = 0          ! record number of slow time step run.
  type(domain2d), pointer :: domain    => NULL()
  type(domainUG), pointer :: UG_domain => NULL()
  type(river_type), save :: River

!--- clock id variable
  integer :: slowclock, bndslowclock, physicsclock, diagclock, riverclock

!--- tracer-related constants, types, and data
character(*), parameter :: trtable='/land_mod/river_tracer' ! name of the field manager tracer table
integer, protected, public :: num_species  ! number of river tracers, public for test_river_solo
integer, parameter :: num_phys = 2 ! number of "physical" tracres: currently they are ice and heat content

type tracer_data_type
  character(fm_field_name_len) :: &
      name        = '', & ! name of the tracer
      units       = '', & ! units of the tracer
      flux_units  = '', & ! units of associated flux
      store_units = ''    ! units of associated storage
  character(fm_string_len)     :: longname = '' ! longname of the species
  real :: &
      t_ref  = 298.0, &
      vf_ref = 0.0,   &
      q10    = 1.0,   &
      kinv   = 1.0
end type

type(tracer_data_type), allocatable :: trdata(:) ! common tracer data

contains


!#####################################################################
  subroutine river_init( land_lon, land_lat, time, dt_fast, land_domain, land_UG_domain, &
                         land_frac, discharge_tol_in, clw_in, csw_in )
    real,            intent(in) :: land_lon(:,:)     ! geographical longitude of cell center
    real,            intent(in) :: land_lat(:,:)     ! geographical latitude of cell center
    type(time_type), intent(in) :: time              ! current time
    type(time_type), intent(in) :: dt_fast           ! fast time step
    type(domain2d),  intent(in), target :: land_domain       ! land domain
    type(domainUG), intent(in), target :: land_UG_domain
    real,            intent(in) :: land_frac(:,:)       ! land area fraction from land model on UG_domain
    real,            intent(in) :: discharge_tol_in, clw_in, csw_in

    integer              :: unit, io_status, ierr, id_restart
    integer              :: sec, day, i, j, i_species
    integer              :: nxc, nyc
    character(len=128)   :: filename
    integer              :: id_lon, id_lat, id_lonb, id_latb
    type(Leo_Mad_trios)   :: DHG_exp            ! downstream equation exponents
    type(Leo_Mad_trios)   :: DHG_coef           ! downstream equation coefficients
    type(Leo_Mad_trios)   :: AAS_exp            ! at-a-station equation exponents

    type(restart_file_type) :: river_restart

    riverclock = mpp_clock_id('update_river'           , CLOCK_FLAG_DEFAULT, CLOCK_SUBCOMPONENT)
    slowclock = mpp_clock_id('update_river_slow'       , CLOCK_FLAG_DEFAULT, CLOCK_ROUTINE)
    bndslowclock = mpp_clock_id('update_river_bnd_slow', CLOCK_FLAG_DEFAULT, CLOCK_ROUTINE)
    physicsclock = mpp_clock_id('river phys'           , CLOCK_FLAG_DEFAULT, CLOCK_ROUTINE)
    diagclock    = mpp_clock_id('river diag'           , CLOCK_FLAG_DEFAULT, CLOCK_ROUTINE)

!--- read namelist -------------------------------------------------
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=river_nml, iostat=io_status)
     ierr = check_nml_error(io_status, 'river_nml')
#else
    if (file_exist('input.nml')) then
      unit = open_namelist_file()
      ierr = 1;
      do while (ierr /= 0)
         read  (unit, nml=river_nml, iostat=io_status, end=10)
         ierr = check_nml_error(io_status,'river_nml')
      enddo
10    continue
      call close_file (unit)
    endif
#endif

!--- write version and namelist info to logfile --------------------
    call log_version(version, module_name, &
    __FILE__)
    unit=stdlog()
    write(unit, river_nml)

    if(.not.do_rivers) return ! do nothing further if the rivers are turned off

!--- check name list variables

    if(diag_freq .le. 0) call mpp_error(FATAL,'river_mod: diag_freq should be a positive integer')

! set up time-related values
    River%time = time
    call get_time(dt_fast, sec, day)
    River%dt_fast = day*sec_in_day+sec

    River%dt_slow = dt_slow
    River%channel_tau = channel_tau

    num_fast_calls = River%dt_slow/River%dt_fast
    River%num_species = num_species
    River%num_c = num_species-num_phys
    River%num_phys = num_phys
    River%i_age = river_tracer_index('age')
    call mpp_error(NOTE,'river_mod: tracer numbers: num_phys='//string(num_phys)//' num_species='//string(num_species))

    if(River%dt_slow .lt. River%dt_fast) call mpp_error(FATAL, &
         'river_mod: river slow time step dt_slow should be no less than land model fast time step dt_fast')

    if ( mod(River%dt_slow,River%dt_fast) .ne. 0  ) call mpp_error(FATAL, &
         'river_mod: river slow time step dt_slow should be multiple of land model fast time step dt_fast')

    discharge_tol = discharge_tol_in
    clw = clw_in
    csw = csw_in

!--- get the domain decompsition, river and land will be on the same grid and have the same domain decomposition.
    domain => land_domain
    UG_domain => land_UG_domain
    call mpp_get_global_domain (domain, xsize=River%nlon, ysize=River%nlat)
    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain   (domain, isd, ied, jsd, jed)
    call mpp_get_UG_compute_domain(UG_domain, lsc, lec)

!---- make sure the halo size is 1
    if( ied-iec .NE. 1 .OR. isc-isd .NE. 1 .OR. jed-jec .NE. 1 .OR. jsc-jsd .NE. 1 ) &
      call mpp_error(FATAL, "river_mod: halo size in four direction should all be 1")

    nxc = iec - isc + 1; nyc = jec - jsc + 1
    !--- make sure land_lon, land_lat is on the compute domain
    if(size(land_lon,1) .NE. nxc .OR. size(land_lon,2) .NE. nyc ) call mpp_error(FATAL, &
        "river_mod: land_lon should be on the compute domain")
    if(size(land_lat,1) .NE. nxc .OR. size(land_lat,2) .NE. nyc ) call mpp_error(FATAL, &
        "river_mod: land_lat should be on the compute domain")

    allocate(missing_rivers        (isc:iec,jsc:jec            ))
    allocate(discharge2ocean_next  (isc:iec,jsc:jec            ))
    allocate(discharge2ocean_next_c(isc:iec,jsc:jec,num_species))
    allocate(id_infloc (0:num_species), id_storage(0:num_species))
    allocate(id_inflow (0:num_species), id_outflow(0:num_species))
    allocate(id_dis    (0:num_species), id_lake_outflow (0:num_species))
    allocate(id_removal(0:num_species), id_stordis(0:num_species), id_run_stor(0:num_species))
    discharge2ocean_next = 0
    discharge2ocean_next_c = 0
    ! IDs of diag fields normalized per cell area
    allocate(id_infloc_c (0:num_species), id_storage_c(0:num_species))
    allocate(id_inflow_c (0:num_species), id_outflow_c(0:num_species))
    allocate(id_dis_c    (0:num_species), id_lake_outflow_c (0:num_species))
    allocate(id_removal_c(0:num_species), id_stordis_c(0:num_species), id_run_stor_c(0:num_species))

!--- read the data from the file river_src_file -- has all static river network data
    call get_river_data(land_lon, land_lat, land_frac)

    missing_rivers = (lnd_sg%landfrac .gt. 0 .and. .not. River%mask)

    i_river_ice  = river_tracer_index('ice')
    i_river_heat = river_tracer_index('het')
    i_river_DOC  = river_tracer_index('doc')
    if (i_river_ice  == NO_TRACER) call mpp_error(FATAL, 'river_mod: required river tracer for ice not found')
    if (i_river_heat == NO_TRACER) call mpp_error(FATAL, 'river_mod: required river tracer for heat not found')

! TODO: get rid of the parameters below in "River" data structure and use trdata
! directly
    River%t_ref  = trdata(num_phys+1:num_species)%t_ref
    River%vf_ref = trdata(num_phys+1:num_species)%vf_ref
    River%q10    = trdata(num_phys+1:num_species)%q10
    River%kinv   = trdata(num_phys+1:num_species)%kinv

!    if (do_age) flux_units(3)                   = 'kg/m2  '
!    if (do_age) store_units(3)                  = 'kg-s/m2'

!--- register diag field
    if(mpp_get_ntile_count(domain)==1) then
       ! grid has just one tile, so we assume that the grid is regular lat-lon
       ! define longitude axes and its edges
       id_lonb = diag_axis_init ( &
            'lonb', lnd_sg%coord_glonb, 'degrees_E', 'X', 'longitude edges', &
            set_name='river', domain2=domain )
       id_lon  = diag_axis_init (                                                &
            'lon',  lnd_sg%coord_glon, 'degrees_E', 'X',  &
            'longitude', set_name='river',  edges=id_lonb, domain2=domain )

       ! define latitude axes and its edges
       id_latb = diag_axis_init ( &
            'latb', lnd_sg%coord_glatb, 'degrees_N', 'Y', 'latitude edges',  &
            set_name='river',  domain2=domain   )
       id_lat = diag_axis_init (                                                &
            'lat',  lnd_sg%coord_glat, 'degrees_N', 'Y', &
            'latitude', set_name='river', edges=id_latb, domain2=domain   )
    else
       id_lon = diag_axis_init ( 'grid_xt', (/(real(i),i=1,River%nlon)/), 'degrees_E', 'X', &
            'T-cell longitude', set_name='river',  domain2=domain, aux='geolon_t' )
       id_lat = diag_axis_init ( 'grid_yt', (/(real(i),i=1,River%nlat)/), 'degrees_N', 'Y', &
            'T-cell latitude', set_name='river',  domain2=domain, aux='geolat_t' )
    endif

    call river_diag_init (id_lon, id_lat)

!--- read restart file
    call get_instance_filename('INPUT/river.res.nc', filename)
    call get_mosaic_tile_file(trim(filename), filename, .false., domain)

    if(file_exist(trim(filename),domain) ) then
        call mpp_error(NOTE, 'river_init : Read restart files '//trim(filename))
        if(new_land_io) then
           id_restart = register_restart_field(river_restart,'river.res.nc','storage', River%storage, domain)
           id_restart = register_restart_field(river_restart,'river.res.nc','discharge2ocean', discharge2ocean_next, domain)
           if (field_exist(filename,'discharge2ocean_c',domain)) then
              id_restart = register_restart_field(river_restart,'river.res.nc','discharge2ocean_c',discharge2ocean_next_c, domain)
              id_restart = register_restart_field(river_restart,'river.res.nc','storage_c', River%storage_c, domain)
           else
              do i_species = 1, num_species
                if (field_exist(filename,'disch2ocn_'//trdata(i_species)%name,domain)) then
                  id_restart = register_restart_field(river_restart,'river.res.nc','disch2ocn_'//trdata(i_species)%name, &
                               discharge2ocean_next_c(:,:,i_species),domain)
                else
                  call mpp_error(WARNING, 'river_init: disch2ocn_'//trim(trdata(i_species)%name)//' does not exist in '//trim(filename))
                endif
                if (field_exist(filename,'storage_'//trdata(i_species)%name,domain)) then
                  id_restart = register_restart_field(river_restart,'river.res.nc','storage_'//trdata(i_species)%name, &
                               River%storage_c(:,:,i_species),domain)
                else
                  call mpp_error(WARNING, 'river_init: storage_'//trim(trdata(i_species)%name)//' does not exist in '//trim(filename))
                endif
              enddo
           endif
           id_restart = register_restart_field(river_restart,'river.res.nc','Omean', River%outflowmean, domain)
           if (field_exist(filename,'depth',domain)) then
              id_restart = register_restart_field(river_restart,'river.res.nc','depth', River%depth, domain, mandatory=.false.)
           endif
           call restore_state(river_restart)
           call free_restart_type(river_restart)
        else
           call read_data(filename,'storage',          River%storage,          domain)
           call read_data(filename,'discharge2ocean',  discharge2ocean_next,   domain)
           if (field_exist(filename,'discharge2ocean_c',domain)) then
              call read_data(filename,'discharge2ocean_c',discharge2ocean_next_c, domain)
              call read_data(filename,'storage_c',        River%storage_c,        domain)
           else
              ! NOTE that the base name of the discharge field was deliberately made
              ! different from the name of the older 3D tracer array, to avoid conflicts
              ! with a tracer is called "C"
              do i_species = 1, num_species
                 if (field_exist(filename,'disch2ocn_'//trdata(i_species)%name,domain)) then
                   call read_data(filename,'disch2ocn_'//trdata(i_species)%name,discharge2ocean_next_c(:,:,i_species), domain)
                 else
                    call mpp_error(WARNING, 'river_init: disch2ocn_'//trim(trdata(i_species)%name)//' does not exist in '//trim(filename))
                 endif
                 if (field_exist(filename,'storage_'//trdata(i_species)%name,domain)) then
                   call read_data(filename,'storage_'//trdata(i_species)%name,River%storage_c(:,:,i_species), domain)
                 else
                    call mpp_error(WARNING, 'river_init: storage_'//trim(trdata(i_species)%name)//' does not exist in '//trim(filename))
                 endif
              enddo
           endif
           call read_data(filename,'Omean',            River%outflowmean,      domain)
           if (field_exist(filename,'depth',domain)) then
                ! call mpp_error(WARNING, 'river_init : Reading field "depth" from '//trim(filename))
                call read_data(filename,'depth',       River%depth,            domain)
           else
                ! call mpp_error(WARNING, 'river_init : "depth" is not present in '//trim(filename))
           endif
        endif
    else
        call mpp_error(NOTE, 'river_init : cold start, set data to 0')
        River%storage    = 0.0
        River%storage_c  = 0.0
        discharge2ocean_next   = 0.0
        discharge2ocean_next_c = 0.0
        if(file_exist(river_Omean_file)) then
           call read_data(river_Omean_file, 'Omean', River%outflowmean, domain)
        else
           River%outflowmean = CONST_OMEAN
        end if
    endif
    River%stordis_c = River%dt_slow * discharge2ocean_next_c/DENS_H2O
    River%stordis   = River%dt_slow *(discharge2ocean_next + &
                                      discharge2ocean_next_c(:,:,1))/DENS_H2O
    where(River%outflowmean .le. outflowmean_min) River%outflowmean=outflowmean_min

    maxtravel = maxval(River%travel)
    call mpp_max(maxtravel)

    call river_physics_init(River, domain, id_lon, id_lat)
    call get_Leo_Mad_params(DHG_exp, DHG_coef, AAS_exp)
    River%o_exp  = 1./ (AAS_exp%on_w + AAS_exp%on_d)
    do j = jsc, jec
       do i = isc, iec
          if ( River%reach_length(i,j) > 0.0) then
              River%o_coef(i,j) = River%outflowmean(i,j) / &
                   ((sinuosity*River%reach_length(i,j))*DHG_coef%on_w*DHG_coef%on_d &
                   *(River%outflowmean(i,j)**(DHG_exp%on_w+DHG_exp%on_d)))**River%o_exp
          endif
       enddo
    enddo
    River%d_exp  = AAS_exp%on_d
    River%d_coef = DHG_coef%on_d                        &
         *(River%outflowmean**(DHG_exp%on_d-AAS_exp%on_d))
    River%w_exp  = AAS_exp%on_w
    River%w_coef = DHG_coef%on_w                        &
         *(River%outflowmean**(DHG_exp%on_w-AAS_exp%on_w))

    num_lake_lev = num_l
    module_is_initialized = .TRUE.

  end subroutine river_init

!#####################################################################
! initialize river tracers
subroutine river_tracers_init()

 integer :: i, m, n
 character(fm_field_name_len) :: name ! name of the river tracer
 character(fm_type_name_len)  :: typ  ! type of the river tracer

 ! number of river tracers in the field table (can be 0)
 m = fm_get_length(trtable)

 ! dump river tracer table
 if(.not.fm_dump_list(trtable, recursive=.TRUE.)) &
    call mpp_error(NOTE, 'river_mod: Cannot dump field list "'//trtable//'"')

 ! allocating more space than absolutely necessary, in case water, and "physical
 ! tracers" (ice and heat) are not present in the user-supplied tracer table
 allocate(trdata(0:m+num_phys))
 ! initialize some parameters of the pre-defined species (water and "physical" tracres)
 trdata(0)%name = 'h2o'; trdata(0)%longname = 'h2o mass'
 trdata(0)%units = 'kg'; trdata(0)%flux_units = 'kg/m2/s'; trdata(0)%store_units = 'kg/m2'

 trdata(1)%name = 'ice'; trdata(1)%longname = 'ice mass'
 trdata(1)%units = 'kg/kg'; trdata(1)%flux_units = 'kg/m2/s'; trdata(1)%store_units = 'kg/m2'

 trdata(2)%name = 'het'; trdata(2)%longname = 'sensible heat content'
 trdata(2)%units = 'K'; trdata(2)%flux_units = 'W/m2'; trdata(2)%store_units = 'J/m2'

 ! read generic parameters of the tracers
 do while (fm_loop_over_list(trtable, name, typ, n))
    ! look for the tracer already in the table
    do i = 0,ubound(trdata,1)
       if (trim(trdata(i)%name)==trim(name)) exit ! found existing slot for this tracer
    enddo
    ! if tracer not found, look for an empty slot in the table
    if (i>=ubound(trdata,1)) then
       do i = 0, ubound(trdata,1)
          if (trim(trdata(i)%name)=='') exit ! found an empty slot
       enddo
    endif
    call read_river_tracer_data(name,trdata(i))
 enddo
 ! finally, calculate the actual number of tracers
 do num_species = ubound(trdata,1),0,-1
    if (trdata(num_species)%name/='') exit ! from loop
 enddo

 ! TODO: read specific tracer parameters. Different tracers might have different parameter sets.

 call print_river_tracer_data(stdout())
 call print_river_tracer_data(stdlog())

end subroutine river_tracers_init

!#####################################################################
function num_river_tracers()
   integer num_river_tracers
   num_river_tracers = num_species
end function num_river_tracers

!#####################################################################
function river_tracer_index(name) result(tr)
   integer :: tr
   character(*), intent(in) :: name

   integer :: i

   tr = NO_TRACER
   do i = 1, num_species
      if (name==trdata(i)%name) then
         tr = i;
         exit
      endif
   enddo
end function river_tracer_index

!#####################################################################
! reads the field_table entry for specific tracers and fills in
! generic tracer parameters
subroutine read_river_tracer_data(name,tr)
  character(*), intent(in) :: name
  type(tracer_data_type), intent(inout) :: tr

  ! ---- local vars
  character(fm_path_name_len)  :: listname
  character(fm_path_name_len)  :: current_list

  current_list = fm_get_current_list()
  if (current_list .eq. ' ') call mpp_error(FATAL, 'river_mod: Could not get the current list')
  listname = trtable//'/'//trim(name)
  if (.not.fm_change_list(listname)) call mpp_error(FATAL,'river_mod: Cannot change field manager list to "'//trim(listname)//'"')

  tr%name        = name
  tr%longname    = fm_util_get_string('long_name',   caller='river_mod', default_value=tr%longname,    scalar=.true.)
  tr%units       = fm_util_get_string('units',       caller='river_mod', default_value=tr%units,       scalar=.true.)
  tr%flux_units  = fm_util_get_string('flux_units',  caller='river_mod', default_value=tr%flux_units,  scalar=.true.)
  tr%store_units = fm_util_get_string('store_units', caller='river_mod', default_value=tr%store_units, scalar=.true.)

#define __PARSE__(v) tr%v = fm_util_get_real(#v, caller='river_mod', default_value=tr%v, scalar=.true.)
  __PARSE__(t_ref)
  __PARSE__(vf_ref)
  __PARSE__(q10)
  __PARSE__(kinv)
#undef __PARSE__

  if (.not.fm_change_list(current_list)) call mpp_error(FATAL,'river_mod: Cannot change field manager list to "'//trim(listname)//'"')
end subroutine read_river_tracer_data

!#####################################################################
! prints a table of tracer data to specified output unit
subroutine print_river_tracer_data(unit)
  integer, intent(in) :: unit

  type(table_printer_type) :: table

  call init_with_headers(table, trdata(:)%name)
  call add_row(table, 'longname', trdata(:)%longname)
  call add_row(table, 'units',    trdata(:)%units)
  call add_row(table, 'flux_units', trdata(:)%flux_units)
  call add_row(table, 'store_units', trdata(:)%store_units)
  call add_row(table, 't_ref', trdata(:)%t_ref)
  call add_row(table, 'vf_ref', trdata(:)%vf_ref)
  call add_row(table, 'q10', trdata(:)%q10)
  call add_row(table, 'kinv', trdata(:)%kinv)

  call print(table,unit)
end subroutine print_river_tracer_data

!#####################################################################
  subroutine update_river ( runoff, runoff_c, land2cplr)
    real, dimension(:,:),   intent(in)  :: runoff
    real, dimension(:,:,:), intent(in)  :: runoff_c
    type(land_data_type), intent(inout) :: land2cplr

    ! --- local vars
    real, dimension(size(runoff,1),size(runoff,2)) :: &
        heat_frac_liq,    & ! fraction of runoff heat in liquid
        discharge_l,      & ! discharge of liquid water to ocean
        discharge_sink      ! container to collect small/negative values for later accounting                                              
    real, dimension(size(runoff,1),size(runoff,2),num_species) ::  &
        discharge_c    ! runoff of tracers accumulated over tiles in cell (including ice and heat)

    integer, save :: n = 0  ! fast time step with each slow time step
    integer       :: i_species
    logical       :: used

    call mpp_clock_begin(riverclock)
    if (.not.do_rivers) then
        call mpp_clock_end(riverclock)  !needed for early exit when do_rivers=.false.
        return
    endif

    discharge_l   = discharge2ocean_next
    discharge_c(:,:,1:num_species) = discharge2ocean_next_c(:,:,1:num_species)
! deplete the discharge storage pools
    River%stordis_c = River%stordis_c &
         - River%dt_fast * discharge_c/DENS_H2O
    River%stordis   = River%stordis   &
         - River%dt_fast *(discharge_l + &
                           discharge_c(:,:,1))/DENS_H2O

!  increment time
    River%Time = increment_time(River%Time, River%dt_fast, 0)
    n = n + 1
!--- accumulate runoff ---------------------
    River%run_stor   = River%run_stor   + runoff
    River%run_stor_c = River%run_stor_c + runoff_c

    if(n == num_fast_calls) then
        call mpp_clock_begin(slowclock)
        call update_river_slow(River%run_stor/real(num_fast_calls), &
             River%run_stor_c(:,:,:)/real(num_fast_calls) )
        call mpp_clock_end(slowclock)
        call mpp_clock_begin(bndslowclock)
        call update_river_bnd_slow
        call mpp_clock_end(bndslowclock)
        n = 0
        River%run_stor = 0
        River%run_stor_c = 0
    endif

    discharge_l = discharge_l/lnd_sg%cellarea
    do i_species = 1, num_species
       discharge_c(:,:,i_species) =  discharge_c(:,:,i_species)/lnd_sg%cellarea
    enddo

    ! pass through to ocean the runoff that was not seen by river module because of land_frac diffs.
    ! need to multiply by gfrac to spread over whole cell
    where (missing_rivers) discharge_l = (runoff-runoff_c(:,:,i_river_ice))*lnd_sg%landfrac
    do i_species = 1, num_species
       where (missing_rivers) &
          discharge_c(:,:,i_species) = runoff_c(:,:,i_species)*lnd_sg%landfrac
    enddo

    ! don't send negatives or insignificant values to ocean. put them in the sink instead.
    ! this code does not seem necessary, and default discharge_tol value should be used.
    discharge_sink = 0.0
    where (discharge_l.le.discharge_tol)
       discharge_sink = discharge_sink + discharge_l
       discharge_l    = 0.0
    end where
    where (discharge_c(:,:,i_river_ice).le.discharge_tol)
       discharge_sink     = discharge_sink + discharge_c(:,:,i_river_ice)
       discharge_c(:,:,i_river_ice) = 0.0
    end where

    ! find phase partitioning ratio for discharge sensible heat flux
    where (discharge_l.gt.0. .or. discharge_c(:,:,i_river_ice).gt.0.)
       heat_frac_liq = clw*discharge_l / (clw*discharge_l+csw*discharge_c(:,:,i_river_ice))
    elsewhere
       heat_frac_liq = 1.0
    end where

    ! scale up fluxes sent to ocean to compensate for non-ocean fraction of discharge cell.
    ! split heat into liquid and solid streams
    where (lnd_sg%landfrac.lt.1.)
       land2cplr%discharge           = discharge_l        / (1-lnd_sg%landfrac)
       land2cplr%discharge_snow      = discharge_c(:,:,i_river_ice) / (1-lnd_sg%landfrac)
       land2cplr%discharge_heat      = heat_frac_liq*discharge_c(:,:,i_river_heat) / (1-lnd_sg%landfrac)
       land2cplr%discharge_snow_heat =               discharge_c(:,:,i_river_heat) / (1-lnd_sg%landfrac) &
                                      - land2cplr%discharge_heat
    end where

#ifdef ZMSDEBUG
    do j=lnd_sg%js,lnd_sg%je
       do i=lnd_sg%is,lnd_sg%ie
          call set_current_point(i, j, 1)
          call check_var_range(land2cplr%discharge(i,j), -10., 10., 'gcell DISCHARGE CHECK', &
                               'Liquid Discharge (mm/s)', WARNING)
          call check_var_range(land2cplr%discharge_snow(i,j), -10., 10., 'gcell DISCHARGE CHECK', &
                               'Snow Discharge (mm/s)', WARNING)
          call check_var_range(land2cplr%discharge_heat(i,j), -1000., 1000., 'gcell DISCHARGE CHECK', &
                               'Discharge Heat (W/m^2/s)', WARNING)
          call check_var_range(land2cplr%discharge_snow_heat(i,j), -1000., 1000., 'gcell DISCHARGE CHECK', &
                               'Discharge Snow Heat (W/m^2/s)', WARNING)
       end do
    end do
#endif

    if (id_dis_liq > 0)  used = send_data (id_dis_liq,  discharge_l, lnd%time)
    if (id_dis_ice > 0)  used = send_data (id_dis_ice,  discharge_c(:,:,i_river_ice), lnd%time)
    if (id_dis_heat > 0) used = send_data (id_dis_heat, discharge_c(:,:,i_river_heat), lnd%time)
    if (id_dis_sink > 0) used = send_data (id_dis_sink, discharge_sink, lnd%time)
    if (id_dis_DOC > 0.and.i_river_DOC/=NO_TRACER) &
                         used = send_data (id_dis_DOC,  discharge_c(:,:,i_river_DOC), lnd%time)

    call mpp_clock_end(riverclock)

  end subroutine update_river

!#####################################################################
  subroutine update_river_slow(runoff, runoff_c)
    real, dimension(:,:),   intent(in)  :: runoff
    real, dimension(:,:,:), intent(in)  :: runoff_c

    real, dimension(isd:ied,jsd:jed) :: &
                             lake_sfc_A, lake_sfc_bot, lake_conn
    real, dimension(isd:ied,jsd:jed,num_lake_lev) :: &
                             lake_wl, lake_ws
    real, dimension(isc:iec,jsc:jec) :: &
                             lake_depth_sill, lake_width_sill, lake_backwater, &
                             lake_backwater_1, &
                             lake_whole_area, &
                             rivr_LMASS,       & ! mass of liquid water in rivers in cell
                             rivr_FMASS,       & ! mass of ice in rivers in cell
                             rivr_MELT,        & ! net mass melt in rivers in cell
                             rivr_HEAT           ! sensible heat content of rivers in cell
    real, dimension(isc:iec,jsc:jec,num_lake_lev) :: &
                             lake_T

    real, dimension(lnd%ls:lnd%le) :: &
                             lake_sfc_A_ug, lake_sfc_bot_ug, lake_conn_ug
    real, dimension(lnd%ls:lnd%le,num_lake_lev) :: &
                             lake_wl_ug, lake_ws_ug
    real, dimension(lnd%ls:lnd%le) :: &
                             lake_depth_sill_ug, lake_width_sill_ug, lake_backwater_ug, &
                             lake_backwater_1_ug, &
                             lake_whole_area_ug
    real, dimension(lnd%ls:lnd%le,num_lake_lev) :: &
                             lake_T_ug

    integer                             :: travelnow, lev, l
    type(Leo_Mad_trios)   :: DHG_exp
    type(Leo_Mad_trios)   :: DHG_coef
    type(Leo_Mad_trios)   :: AAS_exp
    integer i,j,k, i_next, j_next, i_species
    type(land_tile_enum_type)     :: ce    ! land tile enumerator
    type(land_tile_type), pointer :: tile  ! pointer to current tile
    logical :: used
    ! variables for data override
    real, dimension(isc:iec,jsc:jec) :: src_conc, src_flux
    logical :: src_flux_overridden, src_conc_overridden

    slow_step = slow_step + 1

    River%infloc   = River%land_area*runoff  /DENS_H2O
    River%infloc_c = 0
    do i_species = 1, num_species
       River%infloc_c(:,:,i_species) = River%land_area*runoff_c(:,:,i_species)/DENS_H2O
       src_conc = 0.0
       src_flux = 0.0
       call data_override('LND','river_src_flux_'//trdata(i_species)%name, src_flux, River%time, override=src_flux_overridden)
       call data_override('LND','river_src_conc_'//trdata(i_species)%name, src_conc, River%time, override=src_conc_overridden)
       if (src_conc_overridden.OR.src_flux_overridden) then
          where (River%land_area.gt.0)  &
               River%infloc_c(:,:,i_species) = River%infloc*src_conc + src_flux
       endif
    enddo

    River%inflow   = 0
    River%inflow_c = 0
    River%lake_outflow   = 0
    River%lake_outflow_c = 0
    River%disw2o = 0.
    River%disc2o = 0.
    River%melt   = 0.
    lake_sfc_A  = 0
    lake_sfc_bot= 0
    lake_T  = 0
    lake_wl = 0
    lake_ws = 0
    lake_depth_sill  = 0
    lake_width_sill  = 0
    lake_whole_area  = 0
    lake_conn   = 0
    lake_backwater = 0
    lake_backwater_1 = 0
    lake_sfc_A_ug  = 0
    lake_sfc_bot_ug= 0
    lake_T_ug  = 0
    lake_wl_ug = 0
    lake_ws_ug = 0
    lake_depth_sill_ug  = 0
    lake_width_sill_ug  = 0
    lake_whole_area_ug  = 0
    lake_conn_ug   = 0
    lake_backwater_ug = 0
    lake_backwater_1_ug = 0

    ce = first_elmt(land_tile_map, ls=lnd%ls)
    do while(loop_over_tiles(ce, tile, l,k))
       if (.not.associated(tile%lake)) cycle
       if (lake_area_bug) then
          lake_sfc_A_ug (l) = tile%frac * lnd%cellarea(l)
       else
          lake_sfc_A_ug (l) = tile%frac * lnd%area(l)
       endif
       do lev = 1, num_lake_lev
         lake_T_ug (l,lev)   = tile%lake%T(lev)
         lake_wl_ug(l,lev)   = tile%lake%wl(lev)
         lake_ws_ug(l,lev)   = tile%lake%ws(lev)
       enddo
       lake_sfc_bot_ug(l) = (sum(tile%lake%wl(:)+tile%lake%ws(:)) &
                               -tile%lake%wl(1)-tile%lake%ws(1) ) &
                                    / DENS_H2O
       lake_depth_sill_ug(l)  = tile%lake%pars%depth_sill
       lake_width_sill_ug(l)  = tile%lake%pars%width_sill
       lake_whole_area_ug(l)  = tile%lake%pars%whole_area
       lake_conn_ug (l)       = tile%lake%pars%connected_to_next
       lake_backwater_ug(l)   = tile%lake%pars%backwater
       lake_backwater_1_ug(l) = tile%lake%pars%backwater_1
    enddo

!z1l: The following might be changed for performance issue. This might be a temporary solution.
    call mpp_pass_UG_to_SG(lnd%domain, lake_sfc_A_ug, lake_sfc_A)
    call mpp_pass_UG_to_SG(lnd%domain, lake_T_ug, lake_T)
    call mpp_pass_UG_to_SG(lnd%domain, lake_wl_ug, lake_wl)
    call mpp_pass_UG_to_SG(lnd%domain, lake_ws_ug, lake_ws)
    call mpp_pass_UG_to_SG(lnd%domain, lake_sfc_bot_ug, lake_sfc_bot)
    call mpp_pass_UG_to_SG(lnd%domain, lake_depth_sill_ug, lake_depth_sill)
    call mpp_pass_UG_to_SG(lnd%domain, lake_width_sill_ug, lake_width_sill)
    call mpp_pass_UG_to_SG(lnd%domain, lake_whole_area_ug, lake_whole_area)
    call mpp_pass_UG_to_SG(lnd%domain, lake_conn_ug, lake_conn)
    call mpp_pass_UG_to_SG(lnd%domain, lake_backwater_ug, lake_backwater)
    call mpp_pass_UG_to_SG(lnd%domain, lake_backwater_1_ug, lake_backwater_1)
    call mpp_update_domains (lake_sfc_A,  domain)
    call mpp_update_domains (lake_sfc_bot,domain)
    call mpp_update_domains (lake_wl, domain)
    call mpp_update_domains (lake_ws, domain)
    call mpp_update_domains (lake_conn,   domain)
    do i=isc,iec
       do j=jsc,jec
          if (River%i_tocell(i,j)/=NO_RIVER_FLAG) then
             i_next = River%i_tocell(i,j)
             j_next = River%j_tocell(i,j)
          else
             ! to avoid indices out of bounds in the lake_sfc_A check
             i_next = i; j_next=j
          endif

          if (lake_backwater(i,j).gt.0.5 .and. lake_sfc_A(i,j).gt.0. .and. &
                                            lake_sfc_A(i_next,j_next).gt.0. ) then
             ! because of river backwater, lake in this cell relaxes toward level of
             ! lake in next cell downstream. (river depth is still simple function
             ! of discharge though.)
             lake_depth_sill(i,j) = lake_sfc_bot(i_next,j_next) &
              +(lake_wl(i_next,j_next,1)+lake_ws(i_next,j_next,1))/DENS_H2O
          elseif (lake_backwater_1(i,j).gt.0.5) then
             ! to determine depth of backwater, lake at coastal cell has base level
             ! set to river depth in same cell
             lake_depth_sill(i,j) = lake_depth_sill(i,j) + River%depth(i,j)
          elseif (lake_conn(i,j).gt.0.5 ) then
             ! for all but furthest dowstream cell of a multi-cell lake,
             ! relax toward level in next cell (same lake) downstream
             if (lake_conn(i_next,j_next).gt.0.5 .or. all_big_outlet_ctn0) then
                  lake_depth_sill(i,j) = lake_sfc_bot(i_next,j_next) &
                   +(lake_wl(i_next,j_next,1)+lake_ws(i_next,j_next,1))/DENS_H2O
             endif
          elseif (river_impedes_lake) then
              if (lake_width_sill(i,j).lt.0..or.river_impedes_large_lake) then
                  ! lake level in cell relaxes toward river level in cell
                  lake_depth_sill(i,j) = lake_depth_sill(i,j) + River%depth(i,j)
              endif
          endif
       enddo
    enddo

! leftovers from horizontal mixing option, now gone
!call mpp_update_domains (lake_T,  domain)
!call mpp_update_domains (lake_depth_sill, domain)
!call mpp_update_domains (lake_tau, domain)

    travelnow = maxtravel
    do travelnow = maxtravel, 0, -1
       call mpp_clock_begin(physicsclock)
!***************************************************************
       call river_physics_step (River, travelnow, &
         lake_sfc_A, lake_sfc_bot, lake_depth_sill, &
         lake_width_sill, lake_whole_area,         &
         lake_T, lake_wl, lake_ws )
!***************************************************************
       call mpp_clock_end(physicsclock)
    enddo

    call mpp_pass_SG_to_UG(lnd%domain, lake_T, lake_T_ug)
    call mpp_pass_SG_to_UG(lnd%domain, lake_wl, lake_wl_ug)
    call mpp_pass_SG_to_UG(lnd%domain, lake_ws, lake_ws_ug)
   
    ce = first_elmt(land_tile_map, ls=lnd%ls)
    do while(loop_over_tiles(ce, tile, l,k))
       if (.not.associated(tile%lake)) cycle
       do lev = 1, num_lake_lev
         tile%lake%T(lev)  = lake_T_ug (l,lev)
         tile%lake%wl(lev) = lake_wl_ug(l,lev)
         tile%lake%ws(lev) = lake_ws_ug(l,lev)
       enddo
    enddo

    River%outflowmean = River%outflowmean + &
       (River%outflow-River%outflowmean)*River%dt_slow/River%channel_tau
    where(River%outflowmean .le. outflowmean_min) River%outflowmean=outflowmean_min
    call get_Leo_Mad_params(DHG_exp, DHG_coef, AAS_exp)
    do j = jsc, jec
       do i = isc, iec
          if ( River%reach_length(i,j) > 0.0) then
              River%o_coef(i,j) = River%outflowmean(i,j) / &
                   ((sinuosity*River%reach_length(i,j))*DHG_coef%on_w*DHG_coef%on_d &
                   *(River%outflowmean(i,j)**(DHG_exp%on_w+DHG_exp%on_d)))**River%o_exp
          endif
       enddo
    enddo
    River%d_coef = DHG_coef%on_d                        &
         *(River%outflowmean**(DHG_exp%on_d-AAS_exp%on_d))
    River%w_coef = DHG_coef%on_w                        &
         *(River%outflowmean**(DHG_exp%on_w-AAS_exp%on_w))

    River%stordis = River%dt_slow*River%disw2o
    do i_species = 1, num_species
       River%stordis_c(:,:,i_species) = River%dt_slow*River%disc2o(:,:,i_species)
    enddo

    rivr_FMASS = DENS_H2O * (River%storage_c(:,:,1) + River%stordis_c(:,:,1))
    rivr_LMASS = DENS_H2O * (River%storage + River%stordis) - rivr_FMASS
    rivr_MELT  = DENS_H2O *  River%melt / River%dt_slow
    rivr_HEAT  = DENS_H2O * (River%storage_c(:,:,2) + River%stordis_c(:,:,2)) &
                      - hlf*rivr_FMASS

    call mpp_clock_begin(diagclock)
  ! convert area-integrated river outputs to unit-area quantities,
  ! using land area for stores, cell area for fluxes to ocean
    if (id_LWSr > 0) then
       where (lnd_sg%area > 0) &
            rivr_LMASS = rivr_LMASS / lnd_sg%area
       used = send_data (id_LWSr, rivr_LMASS, River%time, mask=lnd_sg%area>0)
    endif
    if (id_FWSr > 0) then
       where (lnd_sg%area > 0) &
            rivr_FMASS = rivr_FMASS / lnd_sg%area
       used = send_data (id_FWSr, rivr_FMASS, River%time, mask=lnd_sg%area>0)
    endif
    if (id_HSr > 0) then
       where (lnd_sg%area > 0) &
            rivr_HEAT = rivr_HEAT / lnd_sg%area
       used = send_data (id_HSr, rivr_HEAT, River%time, mask=lnd_sg%area>0)
    endif
    if (id_meltr > 0) then
       where (lnd_sg%area > 0) &
            rivr_MELT = rivr_MELT / lnd_sg%area
       used = send_data (id_meltr, rivr_MELT, River%time, mask=lnd_sg%area>0)
    end if
    if(mod(slow_step, diag_freq) == 0)  call river_diag(lake_depth_sill)
    call mpp_clock_end(diagclock)

  end subroutine update_river_slow

!#####################################################################

  subroutine update_river_bnd_slow
    integer :: i_species
! note that land_area is not the total area of the cell, but just the land area
! within the cell, so it cannot be used to normalize fluxes to all-ocean cells.
! we need a true cell area for normalization, so river will
! just return mass flux per unit time and let land_model divide by area
    discharge2ocean_next = DENS_H2O*(River%disw2o - River%disc2o(:,:,1))

    do i_species = 1, num_species
       discharge2ocean_next_c(:,:,i_species) = DENS_H2O*River%disc2o(:,:,i_species)
    enddo

  end subroutine update_river_bnd_slow

!#####################################################################

  subroutine river_end
    integer :: outunit ! unit number for stdout

    if(.not.do_rivers) return ! do nothing further if rivers are turned off

!--- write out checksum
    outunit=stdout()
    write(outunit,*)"Chksum for storage ==> ", mpp_chksum(River%storage(isc:iec,jsc:jec))
    write(outunit,*)"Chksum for storage_c ==> ", mpp_chksum(River%storage_c(isc:iec,jsc:jec,:))
    write(outunit,*)"Chksum for discharge2ocean_next ==> ", mpp_chksum(discharge2ocean_next(isc:iec,jsc:jec))
    write(outunit,*)"Chksum for discharge2ocean_next_c ==> ", mpp_chksum(discharge2ocean_next_c(isc:iec,jsc:jec,:))

!--- release memory
    deallocate(discharge2ocean_next, discharge2ocean_next_c,&
         River%run_stor, River%run_stor_c)

    deallocate( River%lon, River%lat)
    deallocate(River%land_area ,     River%basinid        )
    deallocate(River%landfrac )
    deallocate(River%tocell )
    deallocate(River%travel )
    deallocate(River%outflow  )
    deallocate(River%inflow  )
    deallocate(River%lake_outflow)
    deallocate(River%lake_outflow_c)
    deallocate(River%storage        )
    deallocate(River%stordis        )
    deallocate(River%melt           )
    deallocate(River%disw2o        )
    deallocate(River%disc2o        )
    deallocate(River%infloc   )
    deallocate(River%reach_length    )
    deallocate(River%mask        )
    deallocate(River%So        )
    deallocate(River%depth     )
    deallocate(River%width     )
    deallocate(River%vel       )
    deallocate(River%infloc_c ,     River%storage_c ,     River%stordis_c    )
    deallocate(River%inflow_c, River%outflow_c )
    deallocate(River%removal_c )
    deallocate(River%vf_ref,River%t_ref,River%q10,River%kinv)
    deallocate(River%d_coef,River%o_coef,River%w_coef)

    deallocate(trdata)
    module_is_initialized = .FALSE.

  end subroutine river_end


!#####################################################################
  !--- write to restart file
  subroutine save_river_restart(timestamp)
    character(*), intent(in) :: timestamp

    character(len=128) :: filename
    integer :: tr, id_restart
    type(restart_file_type) :: river_restart

    if(.not.do_rivers) return ! do nothing further if rivers are turned off
    if(new_land_io) then
       id_restart = register_restart_field(river_restart,trim(timestamp)//'river.res.nc','storage', River%storage, domain)
       id_restart = register_restart_field(river_restart,trim(timestamp)//'river.res.nc','discharge2ocean', discharge2ocean_next, domain)
       do tr = 1, num_species
          id_restart = register_restart_field(river_restart,trim(timestamp)//'river.res.nc','storage_'//trdata(tr)%name, &
                       River%storage_c(:,:,tr),domain)
          id_restart = register_restart_field(river_restart,trim(timestamp)//'river.res.nc','disch2ocn_'//trdata(tr)%name, &
                       discharge2ocean_next_c(:,:,tr),domain)
       enddo
       id_restart = register_restart_field(river_restart,trim(timestamp)//'river.res.nc','Omean', River%outflowmean, domain)
       id_restart = register_restart_field(river_restart,trim(timestamp)//'river.res.nc','depth', River%depth, domain, mandatory=.false.)
       call save_restart(river_restart)
       call free_restart_type(river_restart)
    else

       filename = 'RESTART/'//trim(timestamp)//'river.res.nc'

       call write_data(filename,'storage', River%storage(isc:iec,jsc:jec), domain)
       call write_data(filename,'discharge2ocean', discharge2ocean_next(isc:iec,jsc:jec), domain)

       !--- write out tracer data
       do tr = 1, num_species
           call write_data(filename,'storage_'//trdata(tr)%name,River%storage_c(isc:iec,jsc:jec,tr), domain)
           call write_data(filename,'disch2ocn_'//trdata(tr)%name,discharge2ocean_next_c(isc:iec,jsc:jec,tr), domain)
       end do
       call write_data(filename,'Omean', River%outflowmean, domain)
       call write_data(filename,'depth', River%depth, domain)
    endif

  end subroutine save_river_restart

!#####################################################################
  subroutine get_river_data(land_lon, land_lat, land_frac)
    real,            intent(in) :: land_lon(isc:,jsc:)  ! geographical lontitude of cell center
    real,            intent(in) :: land_lat(isc:,jsc:)  ! geographical lattitude of cell center
    real,            intent(in) :: land_frac(isc:,jsc:) ! land area fraction of land grid.

    integer                           :: ni, nj, i, j, siz(4), ntiles
    real, dimension(:,:), allocatable :: xt, yt, frac, glon, glat, lake_frac
    integer :: nerrors ! number of errors detected during initialization

    ntiles = mpp_get_ntile_count(domain)

    call field_size(river_src_file, 'basin', siz, domain=domain)
    ni = siz(1)
    nj = siz(2)
    if(ni .NE. River%nlon .OR. nj .NE. River%nlat) call mpp_error(FATAL, &
       "river_mod: size mismatch between river grid and land grid")

    allocate(glon(ni,nj), glat(ni, nj))
    allocate(xt(isc:iec, jsc:jec), yt(isc:iec, jsc:jec), frac(isc:iec, jsc:jec) )
    allocate(lake_frac(isc:iec, jsc:jec))

    if (ntiles == 1) then
        call read_data(river_src_file, 'x', glon, no_domain=.true.)
        call read_data(river_src_file, 'y', glat, no_domain=.true.)
      endif
    call read_data(river_src_file, 'x', xt, domain)
    call read_data(river_src_file, 'y', yt, domain)
    call read_data(river_src_file, 'land_frac', frac, domain)
    call read_data(river_src_file, 'lake_frac', lake_frac, domain)
    !--- the following will be changed when the river data sets is finalized.
    xt = land_lon
    yt = land_lat
!--- transform to radians, since land model grid use radians and compare with land grid.

    allocate(River%lon_1d    (1:ni            ) )
    allocate(River%lat_1d    (1:nj            ) )
    allocate(River%lon       (isc:iec, jsc:jec) )
    allocate(River%lat       (isc:iec, jsc:jec) )
    allocate(River%land_area  (isc:iec, jsc:jec) )
    allocate(River%basinid   (isc:iec, jsc:jec) )
    allocate(River%landfrac  (isc:iec, jsc:jec) )
    allocate(River%mask      (isc:iec, jsc:jec) )
    allocate(River%tocell    (isc:iec, jsc:jec) )
    allocate(River%i_tocell  (isc:iec, jsc:jec) )
    allocate(River%j_tocell  (isc:iec, jsc:jec) )
    allocate(River%travel    (isd:ied, jsd:jed) )
    allocate(River%inflow    (isc:iec, jsc:jec) )
    allocate(River%outflow   (isc:iec, jsc:jec) )
    allocate(River%lake_outflow(isc:iec, jsc:jec) )
    allocate(River%storage   (isc:iec, jsc:jec) )
    allocate(River%stordis   (isc:iec, jsc:jec) )
    allocate(River%run_stor  (isc:iec, jsc:jec) )
    allocate(River%melt      (isc:iec, jsc:jec) )
    allocate(River%disw2o    (isc:iec, jsc:jec) )
    allocate(River%infloc    (isc:iec, jsc:jec))
    allocate(River%reach_length(isc:iec, jsc:jec) )
    allocate(River%So        (isc:iec, jsc:jec) )
    allocate(River%depth     (isc:iec, jsc:jec) )
    allocate(River%width    (isc:iec, jsc:jec) )
    allocate(River%vel      (isc:iec, jsc:jec) )
    allocate(River%infloc_c  (isc:iec, jsc:jec, num_species) )
    allocate(River%storage_c (isc:iec, jsc:jec, num_species) )
    allocate(River%stordis_c (isc:iec, jsc:jec, num_species) )
    allocate(River%run_stor_c (isc:iec, jsc:jec, num_species) )
    allocate(River%outflow_c (isc:iec, jsc:jec, num_species) )
    allocate(River%lake_outflow_c (isc:iec, jsc:jec, num_species) )
    allocate(River%removal_c (isc:iec, jsc:jec, num_species) )
    allocate(River%inflow_c  (isc:iec, jsc:jec, num_species) )
    allocate(River%disc2o    (isc:iec, jsc:jec, num_species))
    allocate(River%d_coef    (isc:iec, jsc:jec) )
    allocate(River%o_coef    (isc:iec, jsc:jec) )
    allocate(River%w_coef    (isc:iec, jsc:jec) )
    allocate(River%outflowmean(isc:iec, jsc:jec) )
    allocate(River%t_ref(num_phys+1:num_species),River%vf_ref(num_phys+1:num_species))
    allocate(River%q10  (num_phys+1:num_species),River%kinv  (num_phys+1:num_species))

    if(ntiles == 1) then   ! lat-lon grid, use actual grid location
       River%lon_1d(:)      = glon(:,1)
       River%lat_1d(:)      = glat(1,:)
    else                   ! cubic grid, use index.
       River%lon_1d(:)      = (/ (i, i=1,River%nlon) /)
       River%lat_1d(:)      = (/ (i, i=1,River%nlat) /)
    end if
    deallocate(glon, glat)

    River%lon(:,:)       = land_lon(:,:)
    River%lat(:,:)       = land_lat(:,:)
!!$    River%landfrac(:,:)  = land_frac(:,:)
    River%landfrac(:,:)  = frac(:,:)
    River%infloc    = 0.0
    River%infloc_c  = 0.0
    River%storage   = 0.0
    River%storage_c = 0.0
    River%stordis   = 0.0
    River%run_stor  = 0.0
    River%stordis_c = 0.0
    River%run_stor_c= 0.0
    River%removal_c = 0.0
    River%depth     = 0.
    River%width     = 0.
    River%vel       = 0.
    River%outflow   = 0.
    River%outflow_c = 0.
    River%inflow    = 0.
    River%inflow_c  = 0.
!--- read the data from the source file
    call read_data(river_src_file, 'tocell', River%tocell, domain)

    where (River%tocell(:,:).eq.  4) River%tocell(:,:)=3
    where (River%tocell(:,:).eq.  8) River%tocell(:,:)=4
    where (River%tocell(:,:).eq. 16) River%tocell(:,:)=5
    where (River%tocell(:,:).eq. 32) River%tocell(:,:)=6
    where (River%tocell(:,:).eq. 64) River%tocell(:,:)=7
    where (River%tocell(:,:).eq.128) River%tocell(:,:)=8

    nerrors = 0
    do j = jsc, jec
    do i = isc, iec
!!$          if(abs(xt(i,j) - land_lon(i,j)) > epsln) call mpp_error(FATAL, &
!!$             "get_river_data: longitude mismatch between river grid and land grid")
!!$          if(abs(yt(i,j) - land_lat(i,j)) > epsln) call mpp_error(FATAL, &
!!$             "get_river_data: latitude mismatch between river grid and land grid")
!!$          if(abs(frac(i,j) - land_frac(i,j)) > epsln) call mpp_error(FATAL, &
!!$             "get_river_data: area fraction mismatch between river grid and land grid")

       ! check that river and land masks match
       if ((frac(i,j)>0).neqv.(land_frac(i,j)>0)) then
          call mpp_error(WARNING,'get_river_data: land and river masks do not match at '//&
               trim(coordinates(i,j)))
          nerrors = nerrors+1
       endif

       ! check that the rivers do not discarge in the middle of the continents
       if ((River%tocell(i,j)==0).and.(land_frac(i,j)>1.0-epsln)) then
          print*,i,j,River%tocell(i,j),land_frac(i,j)
          call mpp_error(WARNING, &
               'get_river_data: river discharges into a land point '&
               //trim(coordinates(i,j))//' where there is no ocean')
          nerrors = nerrors+1
       endif
    end do
    end do

    if (nerrors>0.and.stop_on_mask_mismatch) call mpp_error(FATAL,&
        'get_river_data: river/land mask-related mismatch detected during river data initialization')

    call read_data(river_src_file, 'basin', River%basinid, domain)
    where (River%basinid >0)
       River%mask = .true.
    elsewhere
       River%mask = .false.
    endwhere

    River%travel = 0
    call read_data(river_src_file, 'travel', River%travel(isc:iec,jsc:jec), domain)
    call mpp_update_domains(River%travel, domain)
    call read_data(river_src_file, 'celllength', River%reach_length, domain)
    River%reach_length = River%reach_length * River%landfrac * (1.-lake_frac)
    if (land_area_called_cellarea) then
        call read_data(river_src_file, 'cellarea', River%land_area, domain)
      else
        call read_data(river_src_file, 'land_area', River%land_area, domain)
      endif
!    call read_data(river_src_file, 'So', River%So, domain)
    River%So = 0.0
    where (River%So .LT. 0.0) River%So = Somin

    deallocate(lake_frac)

  end subroutine get_river_data

!#####################################################################

  subroutine river_diag_init(id_lon, id_lat)
    integer, intent(in) :: id_lon  ! ID of land longitude (X) diag axis
    integer, intent(in) :: id_lat  ! ID of land latitude (Y) diag axis

    character(len=11)                :: mod_name = 'river'
    real, dimension(isc:iec,jsc:jec) :: tmp, no_riv
    logical                          :: sent
    integer                          :: i
    integer :: id_geolon_t, id_geolat_t, id_area_land, id_cellarea

! static fields
    id_geolon_t = register_static_field ( mod_name, 'geolon_t', (/id_lon,id_lat/), &
         'longitude of grid cell centers', 'degrees_E', missing_value = -1.0e+20 )
    id_geolat_t = register_static_field ( mod_name, 'geolat_t', (/id_lon,id_lat/), &
         'latitude of grid cell centers', 'degrees_N', missing_value = -1.0e+20 )
    id_cellarea = register_static_field ( mod_name, 'cell_area', (/id_lon, id_lat/), &
         'land cell area', 'm2', missing_value = -1.0e+20 )

    id_area_land = register_static_field ( mod_name, 'area_land', (/id_lon, id_lat/), &
         'land area', 'm2', missing_value = -1.0e+20 )

    ! regular diagnostic fields normalized per land area; values outside of land are zeroed
    ! out. This is the traditional way of saving the river diagnostics.
    ! NOTE that for some fields it is problematic, for example when river discharges into
    ! a grid cell where there is no land, and the value is zeroed-out in the output.
    do i = 0, num_species
      id_inflow(i) = register_diag_field ( mod_name, 'rv_i_'//trim(trdata(i)%name),      &
           (/id_lon, id_lat/), River%Time, 'river inflow, '//trim(trdata(i)%longname),   &
           trdata(i)%flux_units, missing_value=missing )
      id_outflow(i) = register_diag_field ( mod_name, 'rv_o_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river outflow, '//trim(trdata(i)%longname),  &
           trdata(i)%flux_units, missing_value=missing )
      id_dis(i)     = register_diag_field ( mod_name, 'rv_d_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'ocean_discharge, '//trim(trdata(i)%longname),&
           trdata(i)%flux_units, missing_value=missing, area=id_area_land )
      call diag_field_add_attribute(id_dis(i),'cell_methods', 'area: mean')
      id_lake_outflow(i) = register_diag_field ( mod_name, 'rv_l_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'lake outflow, '//trim(trdata(i)%longname), &
           trdata(i)%flux_units, missing_value=missing )
      id_infloc(i) = register_diag_field ( mod_name, 'rv_r_'//trim(trdata(i)%name),      &
           (/id_lon, id_lat/), River%Time, 'local runoff, '//trim(trdata(i)%longname),   &
           trdata(i)%flux_units, missing_value=missing )
      id_removal(i) = register_diag_field ( mod_name, 'rv_m_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river removal, '//trim(trdata(i)%longname),  &
           trdata(i)%flux_units, missing_value=missing )
      id_storage(i) = register_diag_field ( mod_name, 'rv_s_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river storage, '//trim(trdata(i)%longname),  &
           trdata(i)%store_units, missing_value=missing, area=id_area_land )
      call diag_field_add_attribute(id_storage(i),'cell_methods', 'area: mean')
      id_stordis(i) = register_diag_field ( mod_name, 'rv_n_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river discharge lag (numerical) storage, '//trim(trdata(i)%longname), &
           trdata(i)%store_units, missing_value=missing )
      id_run_stor(i) = register_diag_field ( mod_name, 'rv_u_'//trim(trdata(i)%name),    &
           (/id_lon, id_lat/), River%Time, 'river runoff lag (numerical) storage, '//trim(trdata(i)%longname), &
           trdata(i)%store_units, missing_value=missing )
    enddo

    ! register fields normalized per cell area
    do i = 0, num_species
      id_inflow_c(i) = register_diag_field ( mod_name, 'rv_i_c_'//trim(trdata(i)%name),    &
           (/id_lon, id_lat/), River%Time, 'river inflow, '//trim(trdata(i)%longname)//', per unit cell area',   &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea)
      call diag_field_add_attribute(id_inflow_c(i),'cell_methods', 'area: mean')

      id_outflow_c(i) = register_diag_field ( mod_name, 'rv_o_c_'//trim(trdata(i)%name),   &
           (/id_lon, id_lat/), River%Time, 'river outflow, '//trim(trdata(i)%longname)//', per unit cell area',  &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_outflow_c(i),'cell_methods', 'area: mean')

      id_dis_c(i)     = register_diag_field ( mod_name, 'rv_d_c_'//trim(trdata(i)%name),   &
           (/id_lon, id_lat/), River%Time, 'ocean_discharge, '//trim(trdata(i)%longname)//', per unit cell area',&
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_dis_c(i),'cell_methods', 'area: mean')

      id_lake_outflow_c(i) = register_diag_field ( mod_name, 'rv_l_c_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'lake outflow, '//trim(trdata(i)%longname)//', per unit cell area', &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_lake_outflow_c(i),'cell_methods', 'area: mean')

      id_infloc_c(i) = register_diag_field ( mod_name, 'rv_r_c_'//trim(trdata(i)%name),    &
           (/id_lon, id_lat/), River%Time, 'local runoff, '//trim(trdata(i)%longname)//', per unit cell area',   &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_infloc_c(i),'cell_methods', 'area: mean')

      id_removal_c(i) = register_diag_field ( mod_name, 'rv_m_c_'//trim(trdata(i)%name),   &
           (/id_lon, id_lat/), River%Time, 'river removal, '//trim(trdata(i)%longname)//', per unit cell area',  &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_removal_c(i),'cell_methods', 'area: mean')

      id_storage_c(i) = register_diag_field ( mod_name, 'rv_s_c_'//trim(trdata(i)%name),   &
           (/id_lon, id_lat/), River%Time, 'river storage, '//trim(trdata(i)%longname)//', per unit cell area',  &
           trdata(i)%store_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_storage_c(i),'cell_methods', 'area: mean')

      id_stordis_c(i) = register_diag_field ( mod_name, 'rv_n_c_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river discharge lag (numerical) storage, '//trim(trdata(i)%longname)//', per unit cell area', &
           trdata(i)%store_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_stordis_c(i),'cell_methods', 'area: mean')

      id_run_stor_c(i) = register_diag_field ( mod_name, 'rv_u_c_'//trim(trdata(i)%name),    &
           (/id_lon, id_lat/), River%Time, 'river runoff lag (numerical) storage, '//trim(trdata(i)%longname)//', per unit cell area', &
           trdata(i)%store_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_run_stor_c(i),'cell_methods', 'area: mean')
    enddo

    id_lake_depth_sill= register_diag_field ( mod_name, 'rv_dsill', (/id_lon, id_lat/), &
         River%Time, 'effective lake sill depth', 'm', missing_value=missing )
    id_outflowmean   = register_diag_field ( mod_name, 'rv_Qavg', (/id_lon, id_lat/), &
         River%Time, 'long-time average vol. flow', 'm3/s', missing_value=missing )
    id_depth     = register_diag_field ( mod_name, 'rv_depth', (/id_lon, id_lat/), &
         River%Time, 'river flow depth', 'm', missing_value=missing )
    id_width     = register_diag_field ( mod_name, 'rv_width', (/id_lon, id_lat/), &
         River%Time, 'river flow width', 'm', missing_value=missing )
    id_vel       = register_diag_field ( mod_name, 'rv_veloc', (/id_lon, id_lat/), &
         River%Time, 'river flow velocity', 'm/s', missing_value=missing )

    id_LWSr   = register_diag_field ( mod_name, 'LWSr', (/id_lon, id_lat/), &
         River%Time, 'river liquid mass storage', 'kg/m2', missing_value=-1.0e+20, area=id_area_land )
    call diag_field_add_attribute(id_LWSr,'cell_methods', 'area: mean')

    id_FWSr   = register_diag_field ( mod_name, 'FWSr', (/id_lon, id_lat/), &
         River%Time, 'river ice mass storage', 'kg/m2', missing_value=-1.0e+20, area=id_area_land )
    call diag_field_add_attribute(id_FWSr,'cell_methods', 'area: mean')

    id_HSr   = register_diag_field ( mod_name, 'HSr', (/id_lon, id_lat/), &
         River%Time, 'river heat storage', 'J/m2', missing_value=-1.0e+20, area=id_area_land )
    call diag_field_add_attribute(id_HSr,'cell_methods', 'area: mean')

    id_meltr   = register_diag_field ( mod_name, 'meltr', (/id_lon, id_lat/), &
         River%Time, 'melt in river system', 'kg/m2/s', missing_value=-1.0e+20, area=id_area_land )
    call diag_field_add_attribute(id_meltr,'cell_methods', 'area: mean')

    id_dis_liq   = register_diag_field ( mod_name, 'dis_liq', (/id_lon, id_lat/), &
         River%time, 'liquid discharge to ocean', 'kg/(m2 s)', missing_value=-1.0e+20, area=id_cellarea )
    call diag_field_add_attribute(id_dis_liq,'cell_methods', 'area: mean')
    id_dis_ice   = register_diag_field ( mod_name, 'dis_ice', (/id_lon, id_lat/), &
         River%time, 'ice discharge to ocean', 'kg/(m2 s)', missing_value=-1.0e+20, area=id_cellarea )
    call diag_field_add_attribute(id_dis_ice,'cell_methods', 'area: mean')
    id_dis_heat   = register_diag_field ( mod_name, 'dis_heat', (/id_lon, id_lat/), &
         River%time, 'heat of mass discharge to ocean', 'W/m2', missing_value=-1.0e+20, area=id_cellarea )
    call diag_field_add_attribute(id_dis_heat,'cell_methods', 'area: mean')
    id_dis_sink   = register_diag_field ( mod_name, 'dis_sink', (/id_lon, id_lat/), &
         River%time, 'burial rate of small/negative discharge', 'kg/(m2 s)', missing_value=-1.0e+20, area=id_cellarea )
    call diag_field_add_attribute(id_dis_sink,'cell_methods', 'area: mean')
    id_dis_DOC    = register_diag_field ( mod_name, 'dis_DOC', (/id_lon, id_lat/), &
         River%time, 'DOC discharge to ocean', 'kgC/m^2/s', missing_value=-1.0e+20 )
    call diag_field_add_attribute(id_dis_sink,'cell_methods', 'area: mean')

    ! static fields

    id_dx = register_static_field ( mod_name, 'rv_length', (/id_lon, id_lat/), &
           'river reach length', 'm', missing_value=missing )
    id_basin = register_static_field ( mod_name, 'rv_basin', (/id_lon, id_lat/), &
           'river basin id', 'none', missing_value=missing )
    id_So = register_static_field ( mod_name, 'So', (/id_lon, id_lat/), &
           'Slope', 'none', missing_value=missing )
    id_travel = register_static_field ( mod_name, 'rv_trav', (/id_lon, id_lat/), &
           'cells left to travel before reaching ocean', 'none', missing_value=missing )
    id_tocell = register_static_field ( mod_name, 'rv_dir', (/id_lon, id_lat/), &
           'outflow direction code', 'none', missing_value=missing )
    id_no_riv = register_static_field ( mod_name, 'no_riv', (/id_lon, id_lat/), &
         'indicator of land without rivers','unitless', missing_value=-1.0 )

    if(id_geolon_t>0) sent=send_data(id_geolon_t, River%lon*180.0/PI, River%time )
    if(id_geolat_t>0) sent=send_data(id_geolat_t, River%lat*180.0/PI, River%time )
    if(id_area_land>0) sent=send_data(id_area_land, lnd_sg%area, River%time )
    if(id_cellarea>0) sent=send_data(id_cellarea, lnd_sg%cellarea, River%time )

    if (id_dx>0) sent=send_data(id_dx, River%reach_length, River%Time, mask=River%mask )
    if (id_basin>0) then
        tmp = River%basinid(isc:iec,jsc:jec)
        sent=send_data(id_basin, tmp, River%Time, mask=River%mask )
      end if
    if (id_So>0) sent=send_data(id_So, River%So, River%Time, mask=River%mask )
    if (id_travel>0) then
        tmp = River%travel(isc:iec,jsc:jec)
        sent=send_data(id_travel, tmp, River%Time, mask=River%mask )
      end if
    if (id_tocell>0) then
        tmp = River%tocell(isc:iec,jsc:jec)
        sent=send_data(id_tocell, tmp, River%Time, mask=River%mask )
      end if

    no_riv = 0.
    where ( lnd_sg%landfrac .gt. 0 .and. .not. River%mask ) no_riv = 1.
    if ( id_no_riv > 0 ) sent = send_data( id_no_riv, no_riv, River%time )


  end subroutine river_diag_init

!#####################################################################

  subroutine river_diag(lake_depth_sill)
    real, dimension(isc:iec,jsc:jec), intent(in) :: lake_depth_sill
    logical :: used   ! logical for send_data
    real diag_factor  (isc:iec,jsc:jec)
    real diag_factor_2(isc:iec,jsc:jec)
    integer :: tr ! iteratior over river tracers

    diag_factor   = DENS_H2O/lnd_sg%cellarea(:,:)
    diag_factor_2 = 1.0/(lnd_sg%cellarea(:,:)*River%dt_slow)

    if (id_inflow_c(0) > 0) used = send_data (id_inflow_c(0), &
            diag_factor*River%inflow(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_outflow_c(0) > 0) used = send_data (id_outflow_c(0), &
            diag_factor*River%outflow(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_storage_c(0) > 0) used = send_data (id_storage_c(0), &
            diag_factor*River%storage(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_stordis_c(0) > 0) used = send_data (id_stordis_c(0), &
            diag_factor*River%stordis(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_run_stor_c(0) > 0) used = send_data (id_run_stor_c(0), &
            River%dt_fast*River%run_stor(isc:iec,jsc:jec)/lnd_sg%cellarea, River%Time, mask=River%mask )
    if (id_infloc_c(0) > 0) used = send_data (id_infloc_c(0), &
            diag_factor*River%infloc(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_dis_c(0) > 0)    used = send_data (id_dis_c(0), &
            diag_factor*River%disw2o(isc:iec,jsc:jec), River%Time)
    if (id_lake_outflow_c(0) > 0) used = send_data (id_lake_outflow_c(0), &
            diag_factor_2*River%lake_outflow(isc:iec,jsc:jec), River%Time, mask=River%mask )

    do tr = 1, num_species
       if (id_outflow_c(tr) > 0) used = send_data (id_outflow_c(tr), &
         diag_factor*River%outflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_lake_outflow_c(tr) > 0) used = send_data (id_lake_outflow_c(tr), &
         diag_factor_2*River%lake_outflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_inflow_c(tr) > 0) used = send_data (id_inflow_c(tr), &
         diag_factor*River%inflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_storage_c(tr) > 0) used = send_data (id_storage_c(tr), &
         diag_factor*River%storage_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_stordis_c(tr) > 0) used = send_data (id_stordis_c(tr), &
         diag_factor*River%stordis_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_run_stor_c(tr) > 0) used = send_data (id_run_stor_c(tr), &
         River%dt_fast*River%run_stor_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_infloc_c(tr) > 0) used = send_data (id_infloc_c(tr), &
         diag_factor*River%infloc_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_removal_c(tr) > 0) used = send_data (id_removal_c(tr), &
         diag_factor*River%removal_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_dis_c(tr) > 0)    used = send_data (id_dis_c(tr), &
         diag_factor*River%disc2o(isc:iec,jsc:jec,tr), River%Time)
    enddo

    ! recalculate area normalization factors and send data normalized per land area.
    diag_factor = 0.
    diag_factor_2 = 0.
    where (River%land_area(isc:iec,jsc:jec).gt.0.) &
                     diag_factor=DENS_H2O/River%land_area(isc:iec,jsc:jec)
    where (River%land_area(isc:iec,jsc:jec).gt.0.) &
                     diag_factor_2=1./(River%land_area(isc:iec,jsc:jec)*River%dt_slow)

    if (id_inflow(0) > 0) used = send_data (id_inflow(0), &
            diag_factor*River%inflow(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_outflow(0) > 0) used = send_data (id_outflow(0), &
            diag_factor*River%outflow(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_storage(0) > 0) used = send_data (id_storage(0), &
            diag_factor*River%storage(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_stordis(0) > 0) used = send_data (id_stordis(0), &
            diag_factor*River%stordis(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_run_stor(0) > 0) used = send_data (id_run_stor(0), &
            River%dt_fast*River%run_stor(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_infloc(0) > 0) used = send_data (id_infloc(0), &
            diag_factor*River%infloc(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_dis(0) > 0)    used = send_data (id_dis(0), &
            diag_factor*River%disw2o(isc:iec,jsc:jec), River%Time)
    if (id_lake_outflow(0) > 0) used = send_data (id_lake_outflow(0), &
            diag_factor_2*River%lake_outflow(isc:iec,jsc:jec), River%Time, mask=River%mask )

    do tr = 1, num_species
       if (id_outflow(tr) > 0) used = send_data (id_outflow(tr), &
         diag_factor*River%outflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_lake_outflow(tr) > 0) used = send_data (id_lake_outflow(tr), &
         diag_factor_2*River%lake_outflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_inflow(tr) > 0) used = send_data (id_inflow(tr), &
         diag_factor*River%inflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_storage(tr) > 0) used = send_data (id_storage(tr), &
         diag_factor*River%storage_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_stordis(tr) > 0) used = send_data (id_stordis(tr), &
         diag_factor*River%stordis_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_run_stor(tr) > 0) used = send_data (id_run_stor(tr), &
         River%dt_fast*River%run_stor_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_infloc(tr) > 0) used = send_data (id_infloc(tr), &
         diag_factor*River%infloc_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_removal(tr) > 0) used = send_data (id_removal(tr), &
         diag_factor*River%removal_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_dis(tr) > 0)    used = send_data (id_dis(tr), &
         diag_factor*River%disc2o(isc:iec,jsc:jec,tr), River%Time)
    enddo

    if (id_lake_depth_sill > 0) used = send_data (id_lake_depth_sill, &
            lake_depth_sill, River%Time, mask=River%mask )
    if (id_outflowmean > 0) used = send_data (id_outflowmean, &
            River%outflowmean(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_width > 0) used = send_data (id_width, &
            River%width(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_depth > 0) used = send_data (id_depth, &
            River%depth(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_vel > 0) used = send_data (id_vel, &
            River%vel(isc:iec,jsc:jec), River%Time, mask=River%mask )

  end subroutine river_diag

!#####################################################################

  subroutine get_Leo_Mad_params(DHG_exp, DHG_coef, AAS_exp)

    type(Leo_Mad_trios), intent(inout) :: DHG_exp  ! Exponents for downstream equations
    type(Leo_Mad_trios), intent(inout) :: DHG_coef ! Coefficients for downstream equations
    type(Leo_Mad_trios), intent(inout) :: AAS_exp  ! Exponents for at-a-station equations

!!! Exponents for the downstream hydraulic geometry equations
    DHG_exp%on_w = ave_DHG_exp(1)
    DHG_exp%on_d = ave_DHG_exp(2)
    DHG_exp%on_V = ave_DHG_exp(3)

!!! Coefficients for the downstream hydraulic geometry equations
    DHG_coef%on_w = ave_DHG_coef(1)
    DHG_coef%on_d = ave_DHG_coef(2)
    DHG_coef%on_V = ave_DHG_coef(3)

!!! Exponents for the at-a-station hydraulic geometry equations
    AAS_exp%on_w = ave_AAS_exp(1)
    AAS_exp%on_d = ave_AAS_exp(2)
    AAS_exp%on_V = ave_AAS_exp(3)

  end subroutine get_Leo_Mad_params

!#####################################################################

subroutine river_stock_pe(index, value)
integer, intent(in)  :: index
real   , intent(out) :: value ! Domain water (Kg) or heat (Joules)

value = 0.0
if (.not.do_rivers) return

select case(index)
case(ISTOCK_WATER)
  value = DENS_H2O*(sum(River%storage)+sum(River%stordis)) &
        + sum(River%run_stor*River%land_area)*River%dt_fast

case(ISTOCK_HEAT)
! heat stock not yet implemented
  value = 0
case default
! Lnd_stock_pe issues a FATAL error message if index is invalid
end select

end subroutine river_stock_pe

!#####################################################################
! returns string indicating the coordiantes of the point i,j
function coordinates(i,j) result(s); character(128) :: s
   integer, intent(in) :: i,j
   s ='('//trim(string(i))//','//trim(string(j))//')'
   if (lnd%nfaces>1) s=trim(s)//' on cubic sphere face '//string(lnd_sg%face)
end function coordinates

end module river_mod


#ifdef test_river_solo

program river_solo
  use mpp_mod,                  only : mpp_error, mpp_pe, mpp_root_pe, mpp_npes, FATAL
  use mpp_mod,                  only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_domains_mod,          only : mpp_define_layout, mpp_define_domains
  use mpp_domains_mod,          only : mpp_get_compute_domain, domain2d, CYCLIC_GLOBAL_DOMAIN
  use mpp_domains_mod,          only : mpp_get_current_ntile, mpp_get_tile_id
  use mpp_io_mod,               only : mpp_open, MPP_RDONLY, MPP_NETCDF, MPP_SINGLE
  use mpp_io_mod,               only : MPP_ASCII, MPP_OVERWR, mpp_close
  use fms_mod,                  only : fms_init, fms_end, stdlog, open_namelist_file
  use fms_mod,                  only : check_nml_error, close_file, file_exist, stdout, read_data
  use fms_io_mod,               only : fms_io_exit
  use time_manager_mod,         only : time_type, increment_time, set_date, increment_date, set_time
  use time_manager_mod,         only : set_calendar_type, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only : operator(/), operator(-), operator( + ), month_name, get_date
  use diag_manager_mod,         only : diag_manager_init, diag_manager_end
  use river_mod,                only : river_init, river_end, update_river
  use constants_mod,            only : constants_init, PI, radius
  use grid_mod,                 only : get_grid_size, get_grid_cell_vertices
  use grid_mod,                 only : get_grid_cell_centers, get_grid_cell_area, get_grid_comp_area
  use grid_mod,                 only : define_cube_mosaic, get_grid_ntiles
  use river_mod,                only : num_species


  implicit none

  real, parameter       :: CONST_RUNOFF = 200.0

!--- namelist -----------------------------------------

  integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=16)     :: calendar = 'julian'
  integer               :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0
  integer               :: dt_fast     = 0
  integer               ::layout(2) = (/1,0/)
  namelist /river_solo_nml/ current_date, dt_fast, years, months, days, &
       hours, minutes, seconds, calendar, layout

!--------------------------------------------------------------------
  type(time_type)      :: Time, Time_start, Time_end, Run_len, Time_step_fast
  integer              :: nf, num_fast_step, unit
  integer              :: yr,mon,day,hr,min,sec, calendar_type=-1
  integer              :: outunit
  real, allocatable    :: runoff(:,:), discharge(:,:)
  integer              :: initClock, mainClock, termClock, updateClock
!Balaji
  real, allocatable :: runoff_c(:,:,:)
  real, allocatable :: discharge2ocean(:,:)
  real, allocatable :: discharge2ocean_c(:,:,:)
  real, allocatable :: liq(:,:), sol(:,:), mel(:,:), hea(:,:)

  call fms_init

  initClock = mpp_clock_id( 'Initialization' )
  mainClock = mpp_clock_id( 'Main loop' )
  termClock = mpp_clock_id( 'Termination' )
  updateClock = mpp_clock_id( 'update river')

  call mpp_clock_begin(initClock)
  call river_solo_init
  call mpp_clock_end (initClock) !end initialization

  call mpp_clock_begin(mainClock) !begin main loop
  outunit=stdout()
  do nf = 1, num_fast_step
     write(outunit,*)' at river fast time step ', nf
     call mpp_clock_begin(updateClock)
!Balaji
     call update_river( runoff, runoff_c, discharge2ocean, discharge2ocean_c, &
     liq, sol, mel, hea )
!     call update_river ( runoff, discharge )
     call mpp_clock_end (updateClock)
     Time = Time + Time_step_fast
  enddo
  call mpp_clock_end(mainClock)

  call mpp_clock_begin(termClock)
  call river_end
  call diag_manager_end(Time)
  call get_date(Time,yr,mon,day,hr,min,sec)

  if (mpp_pe() == mpp_root_pe()) then
      call mpp_open(unit, 'RESTART/river_solo.res',form=MPP_ASCII,&
           action=MPP_OVERWR,threading=MPP_SINGLE,fileset=MPP_SINGLE,nohdrs=.true.)
      write(unit,*) yr, mon, day, hr, min, sec
      write(unit,*) calendar_type
      call mpp_close(unit)
  endif

  call fms_io_exit
  call fms_end


contains

!#######################################################################
  subroutine river_solo_init

    integer                     :: unit, ierr, io
    integer                     :: ni, nj, npes, isc, iec, jsc, jec, ntiles, tile
    type(domain2d)              :: Domain
    integer, allocatable        :: tile_ids(:)             ! mosaic tile IDs for the current PE
    real,   allocatable         :: lon(:,:), lat(:,:)
    real,   allocatable         :: area_lnd(:,:), area_lnd_cell(:,:), gfrac(:,:)
    integer                     :: date(6)
    character(len=9)            :: month

    call constants_init

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=river_solo_nml, iostat=io)
    ierr = check_nml_error(io, 'river_solo_nml')
#else
    if (file_exist('input.nml')) then
      unit = open_namelist_file ()
      ierr=1
      do while (ierr /= 0)
         read  (unit, nml=river_solo_nml, iostat=io, end=10)
         ierr = check_nml_error (io, 'river_solo_nml')
      enddo
10    continue
      call close_file (unit)
    endif
#endif

    unit=stdlog()
    write(unit, nml= river_solo_nml)

! set the calendar
    if (calendar(1:6) == 'julian') then
        calendar_type = julian
    else if (calendar(1:6) == 'NOLEAP') then
        calendar_type = NOLEAP
    else if (calendar(1:10) == 'thirty_day') then
        calendar_type = THIRTY_DAY_MONTHS
    else if (calendar(1:11) == 'no_calendar') then
        calendar_type = NO_CALENDAR
    else if (calendar(1:1) /= ' ') then
        call mpp_error (FATAL,'==>Error from ocean_solo_mod: invalid namelist value for calendar')
    else
        call mpp_error (FATAL,'==>Error from ocean_solo_mod: no namelist value for calendar')
    endif

! get river_solo restart
    if (file_exist('INPUT/river_solo.res')) then
        call mpp_open(unit,'INPUT/river_solo.res',form=MPP_ASCII,action=MPP_RDONLY)
        read(unit,*) date
        read(unit,*) calendar_type
        call close_file(unit)
    endif

    call set_calendar_type (calendar_type)

    call diag_manager_init

    if (sum(current_date) <= 0) then
        call mpp_error(FATAL,'==>Error from river_solo_mod: no namelist value for current date')
    else
        Time_start  = set_date(current_date(1),current_date(2), current_date(3), &
             current_date(4),current_date(5),current_date(6))
    endif

    if (file_exist('INPUT/river_solo.res')) then
        Time_start =  set_date(date(1),date(2),date(3),date(4),date(5),date(6))
    else
        Time_start = Time_start
        date = current_date
    endif

    Time           = Time_start
    Time_end       = increment_date(Time_start, years, months, days, hours, minutes, seconds)
    Run_len        = Time_end - Time_start
    Time_step_fast = set_time(dt_fast, 0)
    num_fast_step  = Run_len/Time_step_fast

    call mpp_open (unit, 'time_stamp.out', form=MPP_ASCII, action=MPP_OVERWR,threading=MPP_SINGLE)

    month = month_name(current_date(2))
    if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

    call get_date (Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
    month = month_name(date(2))
    if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

    call close_file (unit)

!--- get the land grid and set up domain decomposition
    call get_grid_size('LND', 1, ni, nj)

    npes = mpp_npes()
!--- define domain ------------------------------------------------
    call get_grid_ntiles('LND',ntiles)
    if(layout(1)*layout(2)*ntiles .NE. npes) call mpp_define_layout((/1,ni,1,nj/),npes/ntiles,layout)
    if (ntiles==1) then
       call mpp_define_domains ((/1,ni, 1, nj/), layout, domain, &
            xflags = CYCLIC_GLOBAL_DOMAIN, whalo=1, ehalo=1, shalo=1, nhalo=1, name = 'LAND MODEL')
    else
       call define_cube_mosaic ('LND', domain, layout, halo=1 )
    endif
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)

    allocate(tile_ids(mpp_get_current_ntile(Domain)))
    tile_ids = mpp_get_tile_id(domain)
    tile = tile_ids(1)
    deallocate(tile_ids)

!--- get grid information
    allocate( lon(isc:iec,jsc:jec), lat(isc:iec,jsc:jec) )
    allocate( area_lnd(ni,nj), area_lnd_cell(ni,nj), gfrac(ni,nj) )
    call get_grid_cell_centers ('LND', tile, lon, lat, domain)
!!$    call get_grid_cell_area    ('LND',tile, area_lnd_cell)
!!$    call get_grid_comp_area    ('LND',tile, area_lnd)

    lon = lon * PI/180.
    lat = lat * PI/180.
!!$    gfrac = area_lnd/area_lnd_cell
    gfrac = 1
    npes = mpp_npes()

    call river_init( lon, lat, Time_start, Time_step_fast, Domain, gfrac(isc:iec,jsc:jec)  )

    allocate(runoff_c(isc:iec,jsc:jec,num_species) )
    allocate(runoff(isc:iec,jsc:jec), discharge(isc:iec,jsc:jec) )
    if(file_exist("INPUT/runoff.nc")) then
       call read_data("INPUT/runoff.nc", "runoff", runoff )
    else
       runoff = CONST_RUNOFF
    end if
    if(file_exist("INPUT/runoff.nc")) then
       call read_data("INPUT/runoff.nc", "runoff_c", runoff_c )
    else
       runoff_c = CONST_RUNOFF
    end if
    allocate( discharge2ocean(isc:iec,jsc:jec) )
    allocate( discharge2ocean_c(isc:iec,jsc:jec,num_species) )

  end subroutine river_solo_init

!#####################################################################

end program river_solo

#endif test_river_solo
