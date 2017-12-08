module land_data_mod

use mpp_mod           , only : mpp_get_current_pelist, mpp_pe, mpp_root_pe, &
                               mpp_broadcast, mpp_min, mpp_max, mpp_sum
use constants_mod     , only : PI
use mpp_domains_mod   , only : domain2d, mpp_get_compute_domain, &
     mpp_define_layout, mpp_define_domains, mpp_define_io_domain, &
     mpp_get_current_ntile, mpp_get_tile_id, CYCLIC_GLOBAL_DOMAIN, &
     mpp_get_io_domain, mpp_get_pelist, mpp_get_domain_npes, &
     domainUG, mpp_define_unstruct_domain, mpp_get_UG_domain_tile_id, &
     mpp_get_UG_io_domain, mpp_get_UG_domain_npes, mpp_get_ug_domain_pelist, &
     mpp_get_ug_compute_domain, mpp_get_ug_domain_grid_index, mpp_pass_sg_to_ug, &
     mpp_pass_ug_to_sg, mpp_get_io_domain_UG_layout
use fms_mod           , only : write_version_number, mpp_npes, stdout, &
     file_exist, error_mesg, FATAL, read_data
use fms_io_mod        , only : parse_mask_table
use time_manager_mod  , only : time_type
use grid_mod          , only : get_grid_ntiles, get_grid_size, get_grid_cell_vertices, &
     get_grid_cell_centers, get_grid_cell_area, get_grid_comp_area, &
     define_cube_mosaic
use horiz_interp_mod, only : horiz_interp_type, horiz_interp

implicit none
private

! ==== public interfaces =====================================================
public :: land_data_init
public :: land_data_end
public :: lnd            ! global data
public :: lnd_sg         ! global data related to structured (rectangular) grid
public :: log_version    ! prints version number

public :: atmos_land_boundary_type ! container for information passed from the
                         ! atmosphere to land
public :: land_data_type ! container for information passed from land to
                         ! the atmosphere
public :: land_state_type, land_state_type_sg
public :: horiz_interp_ug
! ==== end of public interfaces ==============================================

! ---- module constants ------------------------------------------------------
character(len=*), parameter :: module_name = 'land_data_mod'
#include "shared/version_variable.inc"

interface prt_max_min_sum
  module procedure prt_max_min_sum_1d
  module procedure prt_max_min_sum_2d
end interface

! ---- types -----------------------------------------------------------------
type :: atmos_land_boundary_type
   ! data passed from the coupler to the surface
   real, dimension(:,:), pointer :: & ! (grid index, tile)
        t_flux    => NULL(), &   ! sensible heat flux, W/m2
        lw_flux   => NULL(), &   ! net longwave radiation flux, W/m2
        lwdn_flux => NULL(), &   ! downward longwave radiation flux, W/m2
        sw_flux   => NULL(), &   ! net shortwave radiation flux, W/m2
        swdn_flux => NULL(), &   ! downward shortwave radiation flux, W/m2
        lprec     => NULL(), &   ! liquid precipitation rate, kg/(m2 s)
        fprec     => NULL(), &   ! frozen precipitation rate, kg/(m2 s)
        tprec     => NULL(), &   ! temperature of precipitation, degK
   ! components of downward shortwave flux, W/m2
        sw_flux_down_vis_dir   => NULL(), & ! visible direct
        sw_flux_down_total_dir => NULL(), & ! total direct
        sw_flux_down_vis_dif   => NULL(), & ! visible diffuse
        sw_flux_down_total_dif => NULL(), & ! total diffuse
   ! derivatives of the fluxes
        dhdt      => NULL(), &   ! sensible w.r.t. surface temperature
        dhdq      => NULL(), &   ! sensible w.r.t. surface humidity
        drdt      => NULL(), &   ! longwave w.r.t. surface radiative temperature
   !
        cd_m      => NULL(), &   ! drag coefficient for momentum, dimensionless
        cd_t      => NULL(), &   ! drag coefficient for tracers, dimensionless
        ustar     => NULL(), &   ! turbulent wind scale, m/s
        bstar     => NULL(), &   ! turbulent buoyancy scale, m/s
        wind      => NULL(), &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot     => NULL(), &   ! height of the bottom atmospheric layer above the surface, m
        drag_q    => NULL(), &   ! product of cd_q by wind
        p_surf    => NULL()      ! surface pressure, Pa

   real, dimension(:,:,:), pointer :: & ! (grid index, tile, tracer)
        tr_flux => NULL(),   &   ! tracer flux, including water vapor flux
        dfdtr   => NULL()        ! derivative of the flux w.r.t. tracer surface value,
                                 ! including evap over surface specific humidity

   integer :: xtype             !REGRID, REDIST or DIRECT
end type atmos_land_boundary_type


type :: land_data_type
   ! data passed from the surface to the coupler
   logical :: pe ! data presence indicator for stock calculations
   real, pointer, dimension(:,:)   :: &  ! (grid index, tile)
        tile_size      => NULL(),  & ! fractional coverage of cell by tile, dimensionless
        t_surf         => NULL(),  & ! ground surface temperature, degK
        t_ca           => NULL(),  & ! canopy air temperature, degK
        albedo         => NULL(),  & ! broadband land albedo [unused?]
        albedo_vis_dir => NULL(),  & ! albedo for direct visible radiation
        albedo_nir_dir => NULL(),  & ! albedo for direct NIR radiation
        albedo_vis_dif => NULL(),  & ! albedo for diffuse visible radiation
        albedo_nir_dif => NULL(),  & ! albedo for diffuse NIR radiation
        rough_mom      => NULL(),  & ! surface roughness length for momentum, m
        rough_heat     => NULL(),  & ! roughness length for tracers and heat, m
        rough_scale    => NULL()     ! topographic scaler for momentum drag, m

   real, pointer, dimension(:,:,:)   :: &  ! (grid index, tile, tracer)
        tr    => NULL()              ! tracers, including canopy air specific humidity

   ! NOTE that in contrast to most of the other fields in this structure, the discharges
   ! hold data per-gridcell, rather than per-tile basis. This, and the order of updates,
   ! have implications for the data reallocation procedure.
   real, pointer, dimension(:,:) :: &  ! (lon, lat)
     discharge           => NULL(),  & ! liquid water flux from land to ocean
     discharge_heat      => NULL(),  & ! sensible heat of discharge (0 C datum)
     discharge_snow      => NULL(),  & ! solid water flux from land to ocean
     discharge_snow_heat => NULL()     ! sensible heat of discharge_snow (0 C datum)

   logical, pointer, dimension(:,:):: &
        mask => NULL()          ! true if land

   integer :: axes(1)           ! ID of diagnostic axes for unstructured grid
   type(domain2D) :: domain     ! structured grid domain
   type(domainUG) :: ug_domain  ! our computation domain
   integer, pointer, dimension(:) :: pelist
end type land_data_type


! land_state_type combines the general information about state of the land model:
! domain, coordinates, time steps, etc. There is only one variable of this type,
! and it is public in this module.
type :: land_state_type_sg
   integer        :: is,ie,js,je ! compute domain boundaries of structured grid domain
   real, allocatable  :: lon (:,:), lat (:,:) ! domain grid center coordinates, radian
   real, allocatable  :: lonb(:,:), latb(:,:) ! domain grid vertices, radian
   real, allocatable  :: area(:,:)  ! land area per grid cell, m2
   real, allocatable  :: cellarea(:,:)  ! grid cell area, m2
   real, allocatable  :: landfrac(:,:)  ! fraction of land in the grid cell
   real, allocatable  :: coord_glon(:), coord_glonb(:) ! longitudes for use in diag axis and such, degrees East
   real, allocatable  :: coord_glat(:), coord_glatb(:) ! latitudes for use in diag axis and such, degrees North

   integer :: face  ! the current mosaic face
   type(domain2D) :: domain ! our domain -- should be the last since it simplifies
                            ! debugging in totalview
end type land_state_type_sg

type :: land_state_type
   integer            :: ls,le       ! domain boundaries of unstructure domain
   integer            :: gs,ge       ! min and max value of grid index ( j*nx+i )
   integer            :: nlon,nlat   ! size of global grid
   type(time_type)    :: dt_fast     ! fast (physical) time step
   type(time_type)    :: dt_slow     ! slow time step

   type(time_type)    :: time        ! current land model time

   real, allocatable  :: area(:)      ! land area per grid cell, m2
   real, allocatable  :: cellarea(:)  ! fraction of land in the grid cell
   real, allocatable  :: landfrac(:)  ! grid cell area, m2
   real, allocatable  :: coord_glon(:), coord_glonb(:) ! longitudes for use in diag axis and such, degrees East
   real, allocatable  :: coord_glat(:), coord_glatb(:) ! latitudes for use in diag axis and such, degrees North
   real, allocatable  :: lon(:), lat(:) ! grid center coordinates, radian
   real, allocatable  :: lonb(:,:), latb(:,:) ! grid vertices, radian

   integer :: nfaces ! number of mosaic faces
   integer :: face  ! the current mosaic face
   integer, allocatable :: pelist(:) ! list of processors that run land model
   integer, allocatable :: io_pelist(:) ! list of processors in our io_domain
   ! if io_domain was not defined, then there is just one element in this
   ! array, and it's equal to current PE
   integer :: io_id     ! suffix in the distributed files.
   logical :: append_io_id ! if FALSE, io_id is not appended to the file names
                           ! (for the case io_layout = 1,1)
   type(domainUG) :: domain ! unstruct grid domain
   integer, allocatable :: i_index(:), j_index(:) ! i,j-index of the unstruture grid on current processor
   integer, allocatable :: l_index(:)             ! l-index (unstructure grid) for idx value.
end type land_state_type


! ---- public module variables -----------------------------------------------
type(land_state_type_sg), protected, save :: lnd_sg
type(land_state_type),    save :: lnd ! It is not protected because lnd%time is updated
                                      ! in other module

! ---- private module variables ----------------------------------------------
logical :: module_is_initialized = .FALSE.


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine log_version(version, modname, filename, tag, unit)
  character(len=*), intent(in) :: version
  character(len=*), intent(in), optional :: &
          modname, filename, tag
  integer, intent(in), optional :: unit
  
  character(512) :: message
  integer :: i
  message=''
  
  if (present(filename)) then
     ! remove the directory part of the name
     i = scan(filename,'/',back=.true.)
     
     message = trim(filename(i+1:))
  endif
  if (present(modname)) then
     message = trim(modname)//' ('//trim(message)//')'
  endif
  call write_version_number (trim(message)//': '//trim(version),tag,unit)
end subroutine log_version

! ============================================================================
subroutine prt_max_min_sum_1d(name, q)
  character(len=*), intent(in)::  name
  real, intent(in)::    q(:)

  real qmin, qmax, qsum

  qmin = minval(q)
  qmax = maxval(q)
  qsum = sum(q)

  call mpp_min(qmin)
  call mpp_max(qmax)
  call mpp_sum(qsum)

  if(mpp_pe() == mpp_root_pe() ) then
        write(*,*) trim(name), ' max = ', qmax, ' min = ', qmin, 'sum = ', qsum
  endif
end subroutine prt_max_min_sum_1d

! ============================================================================
subroutine prt_max_min_sum_2d(name, q)
  character(len=*), intent(in)::  name
  real, intent(in)::    q(:,:)

  real qmin, qmax, qsum

  qmin = minval(q)
  qmax = maxval(q)
  qsum = sum(q)

  call mpp_min(qmin)
  call mpp_max(qmax)
  call mpp_sum(qsum)

  if(mpp_pe() == mpp_root_pe() ) then
        write(*,*) trim(name), ' max = ', qmax, ' min = ', qmin, 'sum = ', qsum
  endif
end subroutine prt_max_min_sum_2d

! ============================================================================
subroutine land_data_init(layout, io_layout, time, dt_fast, dt_slow, mask_table, npes_io_group)
  integer, intent(inout) :: layout(2) ! layout of our domains
  integer, intent(inout) :: io_layout(2) ! layout for land model io
  type(time_type), intent(in) :: &
       time,    & ! current model time
       dt_fast, & ! fast (physical) time step
       dt_slow    ! slow time step
  character(len=*), intent(in) :: mask_table
  integer,          intent(in) :: npes_io_group

  ! ---- local vars
  integer :: nlon, nlat ! size of global grid in lon and lat directions
  integer :: ntiles     ! number of tiles in the mosaic grid
  integer, allocatable :: tile_ids(:) ! mosaic tile IDs for the current PE
  integer :: outunit
  logical :: mask_table_exist
  logical, allocatable :: maskmap(:,:,:)  ! A pointer to an array indicating which
                                          ! logical processors are actually used for
                                          ! the land code. The other logical
                                          ! processors would be all ocean points and
                                          ! are not assigned to actual processors.
                                          ! This need not be assigned if all logical
                                          ! processors are used.

  ! write the version and tag name to the logfile
  call log_version(version, module_name, &
  __FILE__)

  ! define the processor layout information according to the global grid size
  call get_grid_ntiles('LND',ntiles)
  call get_grid_size('LND',1,nlon,nlat)

  mask_table_exist = .false.
  outunit = stdout()
  if(file_exist(mask_table)) then
     mask_table_exist = .true.
     write(outunit, *) '==> NOTE from land_data_init:  reading maskmap information from '//trim(mask_table)
     if(layout(1) == 0 .OR. layout(2) == 0 ) call error_mesg('land_model_init', &
        'land_model_nml layout should be set when file '//trim(mask_table)//' exists', FATAL)

     allocate(maskmap(layout(1), layout(2), ntiles))
     call parse_mask_table(mask_table, maskmap, "Land model")
  endif

  if( layout(1)==0 .AND. layout(2)==0 ) &
       call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes()/ntiles, layout )
  if( layout(1)/=0 .AND. layout(2)==0 )layout(2) = mpp_npes()/(layout(1)*ntiles)
  if( layout(1)==0 .AND. layout(2)/=0 )layout(1) = mpp_npes()/(layout(2)*ntiles)

  if( io_layout(1) == 0 .AND. io_layout(2) == 0 ) io_layout = layout

  ! define land model domain
  if (ntiles==1) then
     if( mask_table_exist ) then
        call mpp_define_domains ((/1,nlon, 1, nlat/), layout, lnd_sg%domain, xhalo=1, yhalo=1,&
             xflags = CYCLIC_GLOBAL_DOMAIN, name = 'LAND MODEL', maskmap=maskmap(:,:,1) )
     else
        call mpp_define_domains ((/1,nlon, 1, nlat/), layout, lnd_sg%domain, xhalo=1, yhalo=1,&
             xflags = CYCLIC_GLOBAL_DOMAIN, name = 'LAND MODEL')
     endif
  else
     if( mask_table_exist ) then
        call define_cube_mosaic ('LND', lnd_sg%domain, layout, halo=1, maskmap=maskmap)
     else
        call define_cube_mosaic ('LND', lnd_sg%domain, layout, halo=1)
     endif
  endif

  ! define io domain
  call mpp_define_io_domain(lnd_sg%domain, io_layout)

  if(mask_table_exist) deallocate(maskmap)

  ! get the domain information of structure domain2d.
  call mpp_get_compute_domain(lnd_sg%domain, lnd_sg%is,lnd_sg%ie,lnd_sg%js,lnd_sg%je)

  ! get the mosaic tile number for this processor: this assumes that there is only one
  ! mosaic tile per PE.
  allocate(tile_ids(mpp_get_current_ntile(lnd_sg%domain)))
  tile_ids = mpp_get_tile_id(lnd_sg%domain)
  lnd_sg%face = tile_ids(1)
  deallocate(tile_ids)

  allocate(lnd_sg%lonb    (lnd_sg%is:lnd_sg%ie+1, lnd_sg%js:lnd_sg%je+1))
  allocate(lnd_sg%latb    (lnd_sg%is:lnd_sg%ie+1, lnd_sg%js:lnd_sg%je+1))
  allocate(lnd_sg%lon     (lnd_sg%is:lnd_sg%ie,   lnd_sg%js:lnd_sg%je))
  allocate(lnd_sg%lat     (lnd_sg%is:lnd_sg%ie,   lnd_sg%js:lnd_sg%je))
  allocate(lnd_sg%area    (lnd_sg%is:lnd_sg%ie,   lnd_sg%js:lnd_sg%je))
  allocate(lnd_sg%cellarea(lnd_sg%is:lnd_sg%ie,   lnd_sg%js:lnd_sg%je))
  allocate(lnd_sg%landfrac(lnd_sg%is:lnd_sg%ie,   lnd_sg%js:lnd_sg%je))
  allocate(lnd_sg%coord_glon(nlon), lnd_sg%coord_glonb(nlon+1))
  allocate(lnd_sg%coord_glat(nlat), lnd_sg%coord_glatb(nlat+1))

  ! initialize coordinates
  call get_grid_cell_vertices('LND',lnd_sg%face,lnd_sg%coord_glonb,lnd_sg%coord_glatb)
  call get_grid_cell_centers ('LND',lnd_sg%face,lnd_sg%coord_glon, lnd_sg%coord_glat)
  call get_grid_cell_area    ('LND',lnd_sg%face,lnd_sg%cellarea, domain=lnd_sg%domain)
  call get_grid_comp_area    ('LND',lnd_sg%face,lnd_sg%area,     domain=lnd_sg%domain)
  lnd_sg%landfrac = lnd_sg%area/lnd_sg%cellarea

  ! set local coordinates arrays -- temporary, till such time as the global arrays
  ! are not necessary
  call get_grid_cell_vertices('LND',lnd_sg%face,lnd_sg%lonb,lnd_sg%latb, domain=lnd_sg%domain)
  call get_grid_cell_centers ('LND',lnd_sg%face,lnd_sg%lon, lnd_sg%lat, domain=lnd_sg%domain)
  ! convert coordinates to radian; note that 1D versions stay in degrees
  lnd_sg%lonb = lnd_sg%lonb*pi/180.0 ; lnd_sg%lon = lnd_sg%lon*pi/180.0
  lnd_sg%latb = lnd_sg%latb*pi/180.0 ; lnd_sg%lat = lnd_sg%lat*pi/180.0

  ! initialize model's time-related parameters
  lnd%time    = time
  lnd%dt_fast = dt_fast
  lnd%dt_slow = dt_slow

  call set_land_state_ug(npes_io_group, ntiles, nlon, nlat)

  module_is_initialized = .TRUE.
end subroutine land_data_init


subroutine set_land_state_ug(npes_io_group, ntiles, nlon, nlat)
  integer,          intent(in) :: npes_io_group, ntiles, nlon, nlat

  ! ---- local vars
  type(domainUG), pointer :: io_domain=>NULL() ! our io_domain
  integer :: n_io_pes ! number of PEs in our io_domain
  integer :: outunit

  !--- variables for unstructure grid domain.
  integer, allocatable :: num_lnd(:), grid_index(:), ntiles_grid(:)
  real,    allocatable :: lnd_area(:,:,:)
  integer              :: i, j, n, l, nland, ug_io_layout

  !-------------------------------------------------------------------
  !   set up for unstructure domain and land state data
  !-------------------------------------------------------------------

  lnd%nlon = nlon; lnd%nlat = nlat
  lnd%nfaces = ntiles
  ! On root pe reading the land_area to decide number of land points.
  allocate(num_lnd(ntiles))

  if(file_exist('INPUT/land_domain.nc', no_domain=.true.)) then
     write(stdout(),*)'set_land_state_ug: reading land information from "INPUT/land_domain.nc" '// &
                      'to use number of land tiles per grid cell for efficient domain decomposition.'
     call read_data('INPUT/land_domain.nc', 'nland_face', num_lnd, no_domain=.true.)
     nland = sum(num_lnd)
     allocate(grid_index(nland))
     allocate(ntiles_grid(nland))
     call read_data('INPUT/land_domain.nc', 'grid_index', grid_index, no_domain=.true.)
     call read_data('INPUT/land_domain.nc', 'grid_ntile', ntiles_grid, no_domain=.true.)
     grid_index = grid_index + 1
  else
     write(stdout(),*)'set_land_state_ug: read land/sea mask from grid file: '// &
                      'number of land tiles per grid cell is not used for domain decomposition'
     if(mpp_pe() == mpp_root_pe()) then
        allocate(lnd_area(nlon,nlat,ntiles))
        do n = 1, ntiles
           call get_grid_comp_area('LND',n,lnd_area(:,:,n))
        enddo
        nland = count(lnd_area(:,:,:) > 0)
        allocate(grid_index(nland))
        l = 0
        num_lnd = 0
        do n = 1,ntiles
           do j = 1, nlat
              do i = 1, nlon
                 if(lnd_area(i,j,n) >0) then
                    l = l+1
                    grid_index(l) = (j-1)*nlon+i
                    num_lnd(n) = num_lnd(n) + 1
                 endif
              enddo
           enddo
        enddo
     endif
     call mpp_broadcast( num_lnd, ntiles, mpp_root_pe() )
     if(mpp_pe() .NE. mpp_root_pe()) then
        nland = sum(num_lnd)
        allocate(grid_index(nland))
     endif 
     call mpp_broadcast( grid_index, nland, mpp_root_pe() )
     allocate(ntiles_grid(nland))
     ntiles_grid = 1
  endif
  call mpp_define_unstruct_domain(lnd%domain, lnd_sg%domain, num_lnd, ntiles_grid, mpp_npes(), &
                                  npes_io_group, grid_index, name = 'LAND MODEL') 
  deallocate(grid_index, num_lnd)
  
  ! set up list of processors  ! set up list of processors for collective io: only the first processor in this
  ! for collective io: only the first processor in this
  ! list actually writes data, the rest just send the data to it.
  io_domain=>mpp_get_UG_io_domain(lnd%domain)
  n_io_pes = mpp_get_UG_domain_npes(io_domain)
  allocate(lnd%io_pelist(n_io_pes))
  call mpp_get_UG_domain_pelist(io_domain,lnd%io_pelist)
  lnd%io_id = mpp_get_UG_domain_tile_id(io_domain)
  ug_io_layout = mpp_get_io_domain_UG_layout(lnd%domain)
  lnd%append_io_id = (ug_io_layout>1)

  ! get the domain information for unstructure domain
  call mpp_get_UG_compute_domain(lnd%domain, lnd%ls,lnd%le)

  !--- get the i,j index of each unstructure grid.
  allocate(grid_index(lnd%ls:lnd%le))
  call mpp_get_UG_domain_grid_index(lnd%domain, grid_index)
  !--- make sure grid_index is monotone increasing.
  do l = lnd%ls+1,lnd%le 
     if(grid_index(l) .LE. grid_index(l-1)) call error_mesg('land_model_init', &
         'grid_index is not monotone increasing', FATAL)
  enddo 
  lnd%gs = grid_index(lnd%ls)
  lnd%ge = grid_index(lnd%le)
  allocate(lnd%l_index(lnd%gs:lnd%ge))
  allocate(lnd%i_index(lnd%ls:lnd%le), lnd%j_index(lnd%ls:lnd%le))
  lnd%l_index = 0
  do l = lnd%ls,lnd%le
     lnd%i_index(l) = mod((grid_index(l)-1), lnd%nlon) + 1
     lnd%j_index(l) = (grid_index(l)-1)/lnd%nlon + 1
     lnd%l_index(grid_index(l)) = l
  enddo
  deallocate(grid_index)
  ! get the mosaic tile number for this processor: this assumes that there is only one
  ! mosaic tile per PE.
  lnd%face = mpp_get_UG_domain_tile_id(lnd%domain)
  allocate(lnd%lonb (lnd%ls:lnd%le, 4))
  allocate(lnd%latb (lnd%ls:lnd%le, 4))
  allocate(lnd%lon  (lnd%ls:lnd%le))
  allocate(lnd%lat  (lnd%ls:lnd%le))
  allocate(lnd%area    (lnd%ls:lnd%le) )
  allocate(lnd%cellarea(lnd%ls:lnd%le) )
  allocate(lnd%landfrac(lnd%ls:lnd%le))
  allocate(lnd%coord_glon(nlon), lnd%coord_glonb(nlon+1))
  allocate(lnd%coord_glat(nlat), lnd%coord_glatb(nlat+1))

  ! initialize coordinates
  call get_grid_cell_vertices('LND',lnd%face,lnd%coord_glonb,lnd%coord_glatb)
  call get_grid_cell_centers ('LND',lnd%face,lnd%coord_glon, lnd%coord_glat)
  call get_grid_cell_area    ('LND',lnd_sg%face,lnd%cellarea, lnd_sg%domain, lnd%domain)
  call get_grid_comp_area    ('LND',lnd_sg%face,lnd%area,     lnd_sg%domain, lnd%domain)
  lnd%landfrac = lnd%area/lnd%cellarea

  ! set local coordinates arrays -- temporary, till such time as the global arrays
  ! are not necessary
  call get_grid_cell_vertices('LND',lnd_sg%face,lnd%lonb,lnd%latb, lnd_sg%domain, lnd%domain)
  call get_grid_cell_centers ('LND',lnd_sg%face,lnd%lon, lnd%lat, lnd_sg%domain, lnd%domain)

!  call mpp_pass_SG_to_UG(lnd%domain, lnd%lon, lnd_sg%lon)
!  call mpp_pass_SG_to_UG(lnd%domain, lnd%lat, lnd_sg%lat)
  ! convert coordinates to radian; note that 1D versions stay in degrees
  lnd%lonb = lnd%lonb*pi/180.0; lnd%lon = lnd%lon*pi/180.0
  lnd%latb = lnd%latb*pi/180.0; lnd%lat = lnd%lat*pi/180.0

  ! initialize the land model processor list
  allocate(lnd%pelist(0:mpp_npes()-1))
  call mpp_get_current_pelist(lnd%pelist)

end subroutine set_land_state_ug


! ============================================================================
subroutine land_data_end()
  module_is_initialized = .FALSE.

end subroutine land_data_end

  subroutine horiz_interp_ug(Interp, data_in, data_out, verbose, &
                             mask_in, mask_out, missing_value, missing_permit, &
                             err_msg, new_missing_handle )
    type (horiz_interp_type),         intent(in) :: Interp
    real, intent(in),             dimension(:,:) :: data_in
    real, intent(out),            dimension(:)   :: data_out
    integer, intent(in),                optional :: verbose
    real, intent(in),   dimension(:,:), optional :: mask_in
    real, intent(out),  dimension(:,:), optional :: mask_out
    real, intent(in),                   optional :: missing_value
    integer, intent(in),                optional :: missing_permit
    character(len=*), intent(out),      optional :: err_msg
    logical, intent(in),                optional :: new_missing_handle

    real, dimension(lnd_sg%is:lnd_sg%ie,lnd_sg%js:lnd_sg%je) :: data_sg

    call horiz_interp(Interp, data_in, data_sg, verbose, &
                                   mask_in, mask_out, missing_value, missing_permit, &
                                   err_msg, new_missing_handle )
    call mpp_pass_sg_to_ug(lnd%domain, data_sg, data_out)                           

  end subroutine horiz_interp_ug

end module land_data_mod
