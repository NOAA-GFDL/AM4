module atmos_cmip_diag_mod

!----------------------------------------------------------------------
!  Module for registering and sending (writing) 3D CMIP diagnostic
!  data on model levels and pressure levels. New vertical axes are 
!  defined from lowest to uppermost level (opposite the model).
!  Prefined pressure axes corresponding to CMIP axes are used.
!  The vertical axis used is specified via the 'module_name' field
!  in the diag_table.
!----------------------------------------------------------------------

use mpp_mod,            only: input_nml_file
use fms_mod,            only: open_namelist_file, check_nml_error, &
                              close_file, stdlog, mpp_pe, mpp_root_pe, &
                              write_version_number, file_exist, &
                              error_mesg, FATAL, WARNING, NOTE, &
                              lowercase, string
use time_manager_mod,   only: time_type
use diag_manager_mod,   only: diag_axis_init, register_diag_field, &
                              send_data, get_diag_field_id, &
                              register_static_field, &
                              diag_axis_add_attribute, &
                              diag_field_add_attribute, &
                              DIAG_FIELD_NOT_FOUND
use diag_data_mod,      only: CMOR_MISSING_VALUE

!----------------------------------------------------------------------

implicit none
private

!----------------------------------------------------------------------

public :: atmos_cmip_diag_init, atmos_cmip_diag_end, &
          register_cmip_diag_field_2d, &
          register_cmip_diag_field_3d, &
          send_cmip_data_3d, &
          query_cmip_diag_id

!----------------------------------------------------------------------

! vertical pressure grids
!    plev   = same as plev17, unless use_extra_levels = .true.
!    plev19 = Table Amon = standard 17 levels + 5, 1 hPa
!    plev8  = daily data
!    plev7h = HighResMIP (6hr time mean, 3hr synoptic)
!    plev3  = used in CMIP5 for 6hrPlev
!    plev23 = plev19 + (7,3,2,0.4 hPa)

real, dimension(23) :: plev23 = &
              (/ 100000., 92500., 85000., 70000., 60000., 50000., &
                  40000., 30000., 25000., 20000., 15000., 10000., &
                   7000.,  5000.,  3000.,  2000.,  1000.,   700., &
                    500.,   300.,   200.,   100.,    40. /)
real, dimension(19) :: plev19 = &
              (/ 100000., 92500., 85000., 70000., 60000., 50000., &
                  40000., 30000., 25000., 20000., 15000., 10000., &
                   7000.,  5000.,  3000.,  2000.,  1000.,   500., &
                    100. /)
real, dimension(8) :: plev8 = &
              (/ 100000., 85000., 70000., 50000., &
                  25000., 10000.,  5000.,  1000. /)
real, dimension(7) :: plev7h = &
              (/  92500., 85000., 70000., 60000., 50000., 25000., 5000. /)
real, dimension(3) :: plev3 = &
              (/  85000., 50000., 25000. /)

!-----------------------------------------------------------------------
!--- namelist ---

logical :: use_extra_levels = .true.  ! use more than the standard
                                      ! 17 pressure levels when possible

logical :: flip_cmip_levels = .true.  ! flip vertical model level output
                                      ! from bottom(surface) to top.

logical :: output_modeling_realm  = .false. ! add modeling_realm attribute
                                            ! to all variables

character(len=64) :: modeling_realm_default = 'atmos' ! default modeling_realm attribute
                                                      ! can be overriden in
                                                      ! register_cmip_diag

integer :: verbose = 1                ! verbose level = 0,1,2

namelist /atmos_cmip_diag_nml/ use_extra_levels, flip_cmip_levels, &
                               output_modeling_realm, modeling_realm_default, &
                               verbose

!-----------------------------------------------------------------------

integer, parameter :: MAXPLEVS = 6  ! max plev sets
integer, dimension(MAXPLEVS) :: num_pres_levs
real,    dimension(MAXPLEVS,50) :: pressure_levels  ! max 50 levels per set

character(len=16) :: mod_name = 'cmip'

! cmip vertical axis names 
!   index -1 = 'levhalf' half model levels
!   index  0 = 'lev'     full model levels
!   index >0 = 'plev*'   pressure levels
character(len=128), dimension(-1:MAXPLEVS) :: cmip_axis_names
integer, dimension(3,-1:MAXPLEVS) :: cmip_axis_data

integer :: area_id

!----------------------------------------------------------------------

!--- store field id for all possible axes
!  index  0 = on model level (either full or half)
!  index >0 = on pressure levels
public :: cmip_diag_id_type
type cmip_diag_id_type
  integer, dimension(0:MAXPLEVS) :: field_id = 0
end type cmip_diag_id_type

!----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

logical :: module_is_initialized=.false.

CONTAINS

!#######################################################################

subroutine atmos_cmip_diag_init ( ak, bk, ptop, axes, Time )
real,    intent(in), dimension(:) :: ak, bk   ! ap,b at model layer interfaces
real,    intent(in)               :: ptop     ! pressure at top level
integer, intent(in)               :: axes(2)  ! x/y axes identifiers
type(time_type), intent(in)       :: Time

!-----------------------------------------------------------------------
! local data

integer :: axes3d(3), k, kk, ind, np, id_plev, num_std_plevs
integer :: nlev
integer :: iunit, ierr, io
logical :: used
character(len=16) :: axis_name
integer :: flip

real                            :: p0
real, dimension(size(ak,1)-1)   :: ap, b, lev
real, dimension(2,size(ak,1)-1) :: ap_bnds, b_bnds, lev_bnds
real, dimension(size(ak,1))   :: levhalf

integer  :: id_lev, id_levhalf, id_nv, id_ap, id_b, &
            id_ap_bnds, id_b_bnds, id_lev_bnds

!-----------------------------------------------------------------------

  if (module_is_initialized) then
    call error_mesg ('atmos_cmip_diag_mod', &
                     'module has already been initialized', WARNING)
    return
  endif

!-----------------------------------------------------------------------
!----- read namelist -----
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=atmos_cmip_diag_nml, iostat=io)
  ierr = check_nml_error (io, 'atmos_cmip_diag_nml')
#else
  if (file_exist('input.nml') ) then
    iunit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read (iunit, nml=atmos_cmip_diag_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'atmos_cmip_diag_nml')
    enddo
10  call close_file (iunit)
  endif
#endif

!----- write version and namelist to log file -----

  iunit = stdlog()
  call write_version_number ( version, tagname )
  if (mpp_pe() == mpp_root_pe()) write (iunit, nml=atmos_cmip_diag_nml)


!-----------------------------------------------------------------------
! axis and area identifiers
  axes3d(1:2) = axes(1:2)
  area_id = get_diag_field_id ('dynamics', 'area')
  if (area_id .eq. DIAG_FIELD_NOT_FOUND) call error_mesg &
        ('atmos_cmip_diag_init', 'diagnostic field "dynamics", '// &
         '"area" is not in the diag_table', NOTE)

!-----------------------------------------------------------------------
! determine the maximum number of standard pressure levels
! first get the pressure (based on ps=1000hPa) at top model level

  if (use_extra_levels) then
    do k = 23, 1, -1
      if (plev23(k) .gt. ptop) then
        num_std_plevs = k
        exit
      endif
    enddo
  else
    num_std_plevs = 17  ! standard ncep 17 levels
  endif

!-----------------------------------------------------------------------
! vertical coordinate variables
! cmip levels are defined from the surface to top (flip model values)

    flip = 1
    if (flip_cmip_levels) flip = -1
    nlev = size(ak,1)-1

    p0 = 100000.
    do k = 1, nlev
      if (flip_cmip_levels) then
        kk = nlev - k + 2
      else
        kk = k
      end if
      ap_bnds(1,k) = ak(kk)
      ap_bnds(2,k) = ak(kk+flip)
      b_bnds(1,k) = bk(kk)
      b_bnds(2,k) = bk(kk+flip)
      ap(k) = (ap_bnds(1,k)+ap_bnds(2,k))*0.5
      b(k) = (b_bnds(1,k)+b_bnds(2,k))*0.5
    enddo
    lev = ap/p0 + b                 ! definition for CMIP purposes
    lev_bnds = ap_bnds/p0 + b_bnds  
    levhalf(1:nlev) = ap_bnds(1,:)/p0 + b_bnds(1,:) ! definition at half levels for CMIP purposes
    levhalf(nlev+1) = ap_bnds(2,nlev)/p0 + b_bnds(2,nlev)

    !---- register new axes ----

    ! at full levels (with bounds attribute)
    id_lev = diag_axis_init('lev', lev, '1.0', 'Z', &
                            'hybrid sigma pressure coordinate', &
                             direction=-1, set_name='cmip', req='lev_bnds')
    call diag_axis_add_attribute (id_lev, 'formula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)')
    call diag_axis_add_attribute (id_lev, 'formula_terms', 'ap: ap b: b ps: ps')
    call diag_axis_add_attribute (id_lev, 'bounds', 'lev_bnds')
    call diag_axis_add_attribute (id_lev, 'standard_name', &
                                  'atmosphere_hybrid_sigma_pressure_coordinate')

    ! at half levels (bounds unknown at top and bottom)
    id_levhalf = diag_axis_init('levhalf', levhalf, '1.0', 'Z', &
                            'hybrid sigma pressure coordinate', &
                             direction=-1, set_name='cmip')
    call diag_axis_add_attribute (id_levhalf, 'standard_name', &
                                  'atmosphere_hybrid_sigma_pressure_coordinate')
    call diag_axis_add_attribute ( id_levhalf, 'formula', 'p(n,k+1/2,j,i) = ap(k+1/2) + b(k+1/2)*ps(n,j,i)')
    call diag_axis_add_attribute ( id_levhalf, 'formula_terms', 'ap: ap_bnds b: b_bnds ps: ps')

    ! vertex number for bounds dimension
    id_nv = diag_axis_init('nv', (/1.,2./), 'none', 'N', 'vertex number', set_name='nv')

    ! register new static variables

    id_ap = register_static_field (mod_name, 'ap', (/id_lev/), &
                                 'vertical coordinate formula term: ap(k)', 'Pa')

    id_b = register_static_field (mod_name, 'b', (/id_lev/), &
                                 'vertical coordinate formula term: b(k)', '1.0')

    id_ap_bnds = register_static_field (mod_name, 'ap_bnds', (/id_nv,id_lev/), &
                                 'vertical coordinate formula term: ap(k+1/2)', 'Pa')

    id_b_bnds = register_static_field (mod_name, 'b_bnds', (/id_nv,id_lev/), &
                                 'vertical coordinate formula term: b(k+1/2)', '1.0')

    id_lev_bnds = register_static_field (mod_name, 'lev_bnds', (/id_nv,id_lev/), &
                                        'hybrid sigma pressure coordinate', '1.0', &
                          standard_name='atmosphere_hybrid_sigma_pressure_coordinate')
    if (id_lev_bnds > 0) then
      call diag_field_add_attribute ( id_lev_bnds, 'formula', 'p(n,k+1/2,j,i) = ap(k+1/2) + b(k+1/2)*ps(n,j,i)')
      call diag_field_add_attribute ( id_lev_bnds, 'formula_terms', 'ap: ap_bnds b: b_bnds ps: ps')
    endif

   ! save static data
   if (id_ap > 0) used = send_data ( id_ap, ap, Time )
   if (id_b  > 0) used = send_data ( id_b , b , Time )
   if (id_ap_bnds  > 0) used = send_data ( id_ap_bnds, ap_bnds, Time )
   if (id_b_bnds   > 0) used = send_data ( id_b_bnds, b_bnds, Time )
   if (id_lev_bnds > 0) used = send_data ( id_lev_bnds, lev_bnds, Time )

   axes3d(3) = id_lev
   cmip_axis_names(0) = 'lev' !mod_name
   cmip_axis_data(:,0) = axes3d

   axes3d(3) = id_levhalf
   cmip_axis_names(-1) = 'levhalf' !mod_name
   cmip_axis_data(:,-1) = axes3d
!-----------------------------------------------------------------------
! loop through all possible pressure level sets
! initialize the pressure axis
! define new 3d grid
! define all 3d state variable on this 3d grid

  do ind = 1, MAXPLEVS
    if (ind .eq. 1) then
      np = num_std_plevs
      pressure_levels(ind,1:np) = plev23(1:np)
      axis_name = 'plev_std'
    else if (ind .eq. 2) then
      np = size(plev19,1)
      pressure_levels(ind,1:np) = plev19
      axis_name = 'plev19'
    else if (ind .eq. 3) then
      np = size(plev8,1)
      pressure_levels(ind,1:np) = plev8
      axis_name = 'plev8'
    else if (ind .eq. 4) then
      np = size(plev3,1)
      pressure_levels(ind,1:np) = plev3
      axis_name = 'plev3'
    else if (ind .eq. 5) then
      np = size(plev7h,1)
      pressure_levels(ind,1:np) = plev7h
      axis_name = 'plev7h'
    else if (ind .eq. 6) then
      np = size(plev23,1)
      pressure_levels(ind,1:np) = plev23
      axis_name = 'plev23'
    endif

    num_pres_levs(ind) = np
    id_plev = diag_axis_init(axis_name, pressure_levels(ind,1:np), &
                             'Pa', 'z', 'pressure', direction=-1, set_name="dynamics")

    axes3d(3) = id_plev
    cmip_axis_names(ind) = trim(axis_name)  !trim(mod_name)//'_'//trim(axis_name)
    cmip_axis_data(:,ind) = axes3d

  enddo

  if (verbose > 0) then
    call error_mesg('atmos_cmip_diag_mod', &
           'cmip_axis_names(-1) = "'//trim(cmip_axis_names(-1))//'"',NOTE)
    do ind = 0, MAXPLEVS
      call error_mesg('atmos_cmip_diag_mod', &
           'cmip_axis_names('//trim(string(ind))//') = "'//trim(cmip_axis_names(ind))//'"',NOTE)
    enddo
  endif

!--- done ---
  module_is_initialized=.true.

!-----------------------------------------------------------------------

end subroutine atmos_cmip_diag_init

!#######################################################################

logical function query_cmip_diag_id (cmip_id, pres)
type(cmip_diag_id_type), intent(in) :: cmip_id
logical, optional,       intent(in) :: pres
integer :: is

  is = 0
  if (present(pres)) then
    if (pres) is = 1
  end if

  query_cmip_diag_id = count(cmip_id%field_id(is:) > 0) .gt. 0

end function query_cmip_diag_id
    
!#######################################################################

integer function register_cmip_diag_field_2d (module_name, field_name, &
                           Time_init, long_name, units, standard_name, &
                           missing_value, interp_method, mask_variant, realm)

  character(len=*), intent(in) :: module_name, field_name
  type(time_type),  intent(in) :: Time_init
  character(len=*), intent(in), optional :: long_name, units, standard_name
  real,             intent(in), optional :: missing_value
  character(len=*), intent(in), optional :: interp_method, realm
  logical         , intent(in), optional :: mask_variant

  real              :: mvalue
  character(len=64) :: modeling_realm
!-----------------------------------------------------------------------
  mvalue = CMOR_MISSING_VALUE; if (present(missing_value)) mvalue = missing_value
  modeling_realm = modeling_realm_default
  if (present(realm)) modeling_realm = realm

  if (output_modeling_realm) then
    register_cmip_diag_field_2d = register_diag_field (module_name, field_name, &
                           cmip_axis_data(1:2,0), Time_init, long_name=long_name, &
                           units=units, standard_name=standard_name, area=area_id, &
                           mask_variant=mask_variant, missing_value=mvalue, &
                           interp_method=interp_method, realm=modeling_realm )
  else
    register_cmip_diag_field_2d = register_diag_field (module_name, field_name, &
                           cmip_axis_data(1:2,0), Time_init, long_name=long_name, &
                           units=units, standard_name=standard_name, area=area_id, &
                           mask_variant=mask_variant, missing_value=mvalue, &
                           interp_method=interp_method)
  endif

!-----------------------------------------------------------------------

end function register_cmip_diag_field_2d

!#######################################################################

function register_cmip_diag_field_3d (module_name, field_name, &
                        Time_init, long_name, units, standard_name, &
                        axis, missing_value, interp_method, mask_variant, &
                        realm)

  character(len=*), intent(in) :: module_name, field_name
  type(time_type),  intent(in) :: Time_init
  character(len=*), intent(in), optional :: long_name, units, standard_name
  real,             intent(in), optional :: missing_value
  character(len=*), intent(in), optional :: axis  ! 'full' or 'half' levels
  character(len=*), intent(in), optional :: interp_method ! for fregrid
  logical         , intent(in), optional :: mask_variant
  character(len=*), intent(in), optional :: realm   ! modeling realm

  type(cmip_diag_id_type) :: register_cmip_diag_field_3d
  integer :: ind, indx, kount
  real    :: mvalue
  character(len=128) :: module_name_table
  character(len=4)   :: vert_axis
  character(len=64)  :: modeling_realm
!-----------------------------------------------------------------------

  mvalue = CMOR_MISSING_VALUE; if (present(missing_value)) mvalue = missing_value
  vert_axis = 'full';          if (present(axis)) vert_axis = lowercase(trim(axis))

  modeling_realm = modeling_realm_default
  if (present(realm)) modeling_realm = realm

  register_cmip_diag_field_3d%field_id = 0

  ! loop thru all axes
  do ind = 0, MAXPLEVS
    indx = ind
    if (ind .eq. 0 .and. vert_axis .eq. 'half') indx = -1

    module_name_table = trim(module_name)
    if (ind .gt. 0) then
      module_name_table = trim(module_name_table)//'_'//trim(cmip_axis_names(ind))
    end if

    ! only register fields that are in the diag_table
    if ( get_diag_field_id(module_name_table, field_name) .ne. DIAG_FIELD_NOT_FOUND ) then

      if (output_modeling_realm) then
        register_cmip_diag_field_3d%field_id(ind) = register_diag_field(module_name_table, field_name,  &
                    cmip_axis_data(:,indx), Time_init, long_name=long_name, units=units, &
                    standard_name=standard_name, area=area_id, mask_variant=mask_variant, &
                    missing_value=mvalue, interp_method=interp_method, realm=modeling_realm)
      else
        register_cmip_diag_field_3d%field_id(ind) = register_diag_field(module_name_table, field_name,  &
                    cmip_axis_data(:,indx), Time_init, long_name=long_name, units=units, &
                    standard_name=standard_name, area=area_id, mask_variant=mask_variant, &
                    missing_value=mvalue, interp_method=interp_method)
      endif

      if (verbose > 0) call error_mesg('atmos_cmip_diag_mod', &
         'register cmip diag: module='//trim(module_name_table)//', field='//trim(field_name)// &
         ', field_id='//trim(string(register_cmip_diag_field_3d%field_id(ind))),NOTE)

    else if (verbose > 1) then
      ! for debugging purposes
      call error_mesg('atmos_cmip_diag_mod','NOT registering cmip diag: module='// &
                         trim(module_name_table)//', field='//trim(field_name),NOTE)
    endif
  enddo

  if (verbose > 1) then
    kount = count(register_cmip_diag_field_3d%field_id > 0)
    if (query_cmip_diag_id(register_cmip_diag_field_3d)) then
      call error_mesg('atmos_cmip_diag_mod','query_cmip_diag_id=TRUE, module='// &
        trim(module_name_table)//', field='//trim(field_name)//', kount='//trim(string(kount)),NOTE)
    else
      call error_mesg('atmos_cmip_diag_mod','query_cmip_diag_id=FALSE, module='// &
        trim(module_name_table)//', field='//trim(field_name)//', kount='//trim(string(kount)),NOTE)
    endif
  endif
    
!-----------------------------------------------------------------------

end function register_cmip_diag_field_3d
 
!#######################################################################

logical function send_cmip_data_3d (cmip_id, field, Time, is_in, js_in, ks_in, phalf, mask, rmask, opt, ext)

  type(cmip_diag_id_type),   intent(in) :: cmip_id
  real, dimension(:,:,:),    intent(in) :: field
  type(time_type),           intent(in), optional :: Time
  integer,                   intent(in), optional :: is_in, js_in, ks_in
  real,    dimension(:,:,:), intent(in), optional :: phalf, rmask
  logical, dimension(:,:,:), intent(in), optional :: mask
  integer,                   intent(in), optional :: opt  ! if opt /= 0 then phalf(i,k,j)
  logical,                   intent(in), optional :: ext

  integer :: ind, id, np, ke
  real, allocatable :: pdat(:,:,:)

!-----------------------------------------------------------------------

  if (.not.module_is_initialized) call error_mesg ('atmos_cmip_diag_mod', &
                               'module has not been initialized', FATAL)

  if (present(ks_in)) then
    if (ks_in .ne. 1) call error_mesg ('atmos_cmip_diag_mod', &
             'subroutine send_cmip_data_3d does not support optional arg "ks_in"', FATAL)
  endif

  if (present(rmask) .and. present(mask)) call error_mesg('atmos_cmip_diag_mod', &
                                   'rmask and mask can not both be present',FATAL)

  send_cmip_data_3d = .false.
  
  ! loop thru all axes
  do ind = 0, MAXPLEVS

    if (cmip_id%field_id(ind) > 0) then
      id = cmip_id%field_id(ind)

    ! pressure level interpolation if "phalf" is present

      ! pressure level interpolation when ind > 0
      if (ind > 0) then
        if (.not.present(phalf)) then
          cycle ! silently skip?
        endif
        if (present(rmask) .or. present(mask)) call error_mesg('atmos_cmip_diag_mod', &
                               'rmask or mask not allowed with pressure interpolation',FATAL)
        np = num_pres_levs(ind)
        allocate(pdat(size(field,1),size(field,2),np))
        call interpolate_vertical_3d (pressure_levels(ind,1:np), phalf, field, pdat, opt=opt, ext=ext)
        send_cmip_data_3d = send_data(id, pdat, Time, is_in=is_in, js_in=js_in, ks_in=ks_in)
        deallocate(pdat)

      else
      ! save data on model levels (flip data)
        if (flip_cmip_levels) then
           ke = size(field,3)
           if (.not.present(mask) .and. .not.present(rmask)) then
             send_cmip_data_3d = send_data(id, field(:,:,ke:1:-1), Time, &
                                      is_in=is_in, js_in=js_in, ks_in=ks_in)
           else if (present(mask) .and. .not.present(rmask)) then
             send_cmip_data_3d = send_data(id, field(:,:,ke:1:-1), Time, &
                                      is_in=is_in, js_in=js_in, ks_in=ks_in, mask=mask(:,:,ke:1:-1))
           else if (.not.present(mask) .and. present(rmask)) then
             send_cmip_data_3d = send_data(id, field(:,:,ke:1:-1), Time, &
                                      is_in=is_in, js_in=js_in, ks_in=ks_in, rmask=rmask(:,:,ke:1:-1))
           endif
         else
           send_cmip_data_3d = send_data(id, field(:,:,:), Time, &
                                      is_in=is_in, js_in=js_in, ks_in=ks_in, mask=mask, rmask=rmask)
         endif
      endif
    else
      send_cmip_data_3d = .false.
    endif
  enddo

!-----------------------------------------------------------------------

end function send_cmip_data_3d

!#######################################################################

subroutine atmos_cmip_diag_end

! do nothing, no way to unregister diag fields

end subroutine atmos_cmip_diag_end

!#######################################################################

subroutine dealloc_cmip_diag_id_type (cmip_id)
class(cmip_diag_id_type), intent(inout) :: cmip_id

 !deallocate(cmip_id%field_id)

end subroutine dealloc_cmip_diag_id_type

!#######################################################################
! wrapper for different vertical interpolation routines
! opt = 0  for standard indexing of peln(i,j,k) -- this is the default
! opt /= 0 for FV-core indexing of peln(i,k,j)
! ext = flag to extrapolate data below (and above) input data (default: false)

subroutine interpolate_vertical_3d (plev, peln, a, ap, opt, ext)
  real,    intent(in),  dimension(:)     :: plev  ! target p-levels
  real,    intent(in),  dimension(:,:,:) :: peln  ! log(phalf), model half levels
  real,    intent(in),  dimension(:,:,:) :: a     ! input data
  real,    intent(out), dimension(:,:,:) :: ap    ! output data on p-levels
  integer, intent(in), optional          :: opt   ! peln indexing
  logical, intent(in), optional          :: ext   ! extrapolate?

  integer :: im, jm, km, kp
  integer :: iopt

  iopt = 0; if (present(opt)) iopt = opt
  im = size(a,1)
  jm = size(a,2)
  km = size(a,3)
  kp = size(plev,1)

  if (iopt .eq. 0) then
    if (size(peln,2).eq.jm .and. size(peln,3).eq.km+1) then
      call interpolate_vertical (im, jm, km, kp, plev, peln, a, ap)  ! peln(im,jm,km+1)
   !else if (size(peln,2).eq.jm .and. size(peln,3).eq.km) then
   !  call interpolate_vertical_half (im, jm, km, kp, plev, peln, a, ap, ext)  ! peln(im,jm,km)
    else
      call error_mesg('atmos_cmip_diag_mod','invalid indexing option and/or array sizes',FATAL)
    endif

  else
    if (size(peln,3).eq.jm .and. size(peln,2).eq.km+1) then
      call interpolate_vertical_fv (im, jm, km, kp, plev, peln, a, ap)  ! peln(im,km+1,jm)
    else if (size(peln,3).eq.jm .and. size(peln,2).eq.km) then
      call interpolate_vertical_half_fv (im, jm, km, kp, plev, peln, a, ap, ext)  ! peln(im,km,jm)
    else
      call error_mesg('atmos_cmip_diag_mod','invalid indexing option and/or array sizes',FATAL)
    endif

  endif

end subroutine interpolate_vertical_3d

!#######################################################################
! a    (im, jm, km  )  <-- input data on FULL model levels
! peln (im, jm, km+1)  <-- standard indexing (i,j,k)
! km = number of FULL levels

subroutine interpolate_vertical (im, jm, km, np, plev, peln, a, ap, ext)
  integer, intent(in)                         :: im, jm, km, np
  real,    intent(in),  dimension(np)         :: plev  ! target p-levels
  real,    intent(in),  dimension(im,jm,km+1) :: peln ! log(phaf), model half levels
  real,    intent(in),  dimension(im,jm,km)   :: a     ! input data on model levels
  real,    intent(out), dimension(im,jm,np)   :: ap    ! output data on p-levels
  logical, intent(in), optional               :: ext

  real    :: pm(km), logp
  integer :: i, j, k, kp
  logical :: extrap

  extrap = .false.; if (present(ext)) extrap = ext

  do kp = 1, np
    logp = log(plev(kp))

    do j = 1, jm
      do i = 1, im
        pm = 0.5*(peln(i,j,1:km)+peln(i,j,2:km+1))
        include "atmos_cmip_interp.inc"
      enddo
    enddo
  enddo

end subroutine interpolate_vertical

!#######################################################################
! a    (im, jm,   km)  <-- input data on FULL model levels
! peln (im, km+1, jm)  <-- FV core indexing
! km = number of FULL levels

subroutine interpolate_vertical_fv (im, jm, km, np, plev, peln, a, ap, ext)
  integer, intent(in)                         :: im, jm, km, np
  real,    intent(in),  dimension(np)         :: plev  ! target p-levels
  real,    intent(in),  dimension(im,km+1,jm) :: peln ! log(phaf), model half levels
  real,    intent(in),  dimension(im,jm,km)   :: a     ! input data on model levels
  real,    intent(out), dimension(im,jm,np)   :: ap    ! output data on p-levels
  logical, intent(in), optional               :: ext

  real    :: pm(km), logp
  integer :: i, j, k, kp
  logical :: extrap

  extrap = .false.; if (present(ext)) extrap = ext

  do kp = 1, np
    logp = log(plev(kp))

    do j = 1, jm
      do i = 1, im
        pm = 0.5*(peln(i,1:km,j)+peln(i,2:km+1,j))
        include "atmos_cmip_interp.inc"
      enddo
    enddo
  enddo

end subroutine interpolate_vertical_fv

!#######################################################################
! a    (im, jm, km)  <-- input data on HALF model levels
! peln (im, km, jm)  <-- FV core indexing
! km = number of HALF levels

subroutine interpolate_vertical_half_fv (im, jm, km, np, plev, peln, a, ap, ext)

  integer, intent(in)                       :: im, jm, km, np
  real,    intent(in),  dimension(np)       :: plev  ! target p-levels
  real,    intent(in),  dimension(im,km,jm) :: peln  ! log(phaf), model half levels
  real,    intent(in),  dimension(im,jm,km) :: a     ! input data on model HALF levels
  real,    intent(out), dimension(im,jm,np) :: ap    ! output data on p-levels
  logical, intent(in), optional             :: ext

  real    :: pm(km), logp
  integer :: i, j, k, kp
  logical :: extrap

  extrap = .false.; if (present(ext)) extrap = ext

  do kp = 1, np
    logp = log(plev(kp))

    do j = 1, jm
      do i = 1, im
        pm = peln(i,:,j)
        include "atmos_cmip_interp.inc"
      enddo
    enddo
  enddo

end subroutine interpolate_vertical_half_fv

!#######################################################################

end module atmos_cmip_diag_mod

