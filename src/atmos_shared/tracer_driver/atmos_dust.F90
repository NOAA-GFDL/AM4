module atmos_dust_mod
! <DESCRIPTION>
!   This module evaluates the change of mass mixing ratio for mineral dust
!   particles due to their emission from preferential sources, and the removal
!   by gravitational settling. The dust particles are transported as dry
!   particles. No hygroscopic growth is considered.
!   The size distribution of sea salt ranges from 0.1 to 10 um (dry radius)
!   and is divided into 5 bins. For each bin, the volume size distribution
!   dV/dlnr is considered constant.
! </DESCRIPTION>
! <CONTACT EMAIL="Paul.Ginoux@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use        constants_mod, only : PI, GRAV, RDGAS, DENS_H2O, PSTD_MKS, WTMAIR
use              mpp_mod, only : input_nml_file 
use              fms_mod, only : write_version_number, mpp_pe,  mpp_root_pe, &
                                 open_namelist_file, close_file, file_exist, &
                                 check_nml_error, error_mesg,  &
                                 stdlog, stdout, string, lowercase, &
                                 NOTE, FATAL
use     time_manager_mod, only : time_type
use     diag_manager_mod, only : send_data, register_diag_field
use  atmos_cmip_diag_mod, only : register_cmip_diag_field_2d
use   tracer_manager_mod, only : get_number_tracers, get_tracer_index, &
                                 get_tracer_names, set_tracer_atts, & 
                                 query_method, NO_TRACER
use    field_manager_mod, only : parse, MODEL_ATMOS, MODEL_LAND
use atmos_tracer_utilities_mod, only : wet_deposition, dry_deposition
use interpolator_mod,      only: interpolate_type, interpolator_init, &
                                 obtain_interpolator_time_slices, &
                                 unset_interpolator_time_flag, &
                                 interpolator, interpolator_end, &
                                 CONSTANT, INTERP_WEIGHTED_P

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_dust_sourcesink, atmos_dust_init, atmos_dust_end, &
        atmos_dust_time_vary, atmos_dust_endts, is_dust_tracer, &
        atmos_dust_flux_init, atmos_dust_gather_data,           &
        dust_has_surf_setl_flux, get_dust_surf_setl_flux,       &
        atmos_dust_wetdep_flux_set,atmos_dust_drydep_flux_set 

public do_dust, n_dust_tracers, dust_tracers

character(len=7), parameter :: module_name = 'tracers'

! data type to hold the individual dust tracer parameters
type :: dust_data_type
   character(32) :: name = '' ! name of the tracer
   integer       :: tr   = NO_TRACER ! index of this dust tracer in the atmos tracer array
   real          :: ra=0.0, rb=0.0   ! boundaries of the size distribution bin, m
   real          :: dustref = 0.0    ! effective radius of the dry dust particles, m
   real          :: dustden = 2650.0 ! density of dry dust particles, kg/m3
   real          :: frac_s = 0.0     ! fraction of total source for this bin
   logical       :: do_surf_exch = .FALSE. ! if true, dust gravitational sedimentation
   ! contributes to the flux between the atmos bottom and the surface in the flux exchange.
   ! Otherwise, it contributes the explicit dust tendency at the bottom layer.
   real, pointer :: dust_setl(:,:) => NULL() ! sedimentation flux at the bottom of the atmos
   real, pointer :: dsetl_dtr(:,:) => NULL() ! derivative of the sedimentation flux w.r.t. dust concentration.
   ! diagnostic IDs
   integer       :: id_dust_emis = -1, id_dust_setl = -1
end type dust_data_type

logical :: do_dust = .FALSE.
integer, save   :: ind_dry_dep_lith_dust_flux = 0
integer, save   :: ind_wet_dep_lith_dust_flux = 0
real, allocatable :: dry_dep_lith_dust_flux(:,:)
real, allocatable :: wet_dep_lith_dust_flux(:,:)

! ---- module data ----
logical :: module_is_initialized = .FALSE.
logical :: do_emission = .TRUE. ! indicates that dust emission is done on atmos. side 
real, save :: u_ts
real, save :: ch
integer :: n_dust_tracers = 0 ! number of dust tracers
type(dust_data_type), allocatable :: dust_tracers(:) ! parameters for specific dust tracers
type(interpolate_type),save       :: dust_source_interp
! ---- identification numbers for diagnostic fields ----
integer :: id_dust_source, id_dust_emis, id_dust_ddep
integer :: id_emidust, id_drydust ! cmip

!---------------------------------------------------------------------
!-------- namelist  ---------
character(len=32)  :: dust_source_filename = 'dust_source_1x1.nc'
character(len=32)  :: dust_source_name(1) = 'source'
real               :: uthresh=-999.
real               :: coef_emis =-999.
logical            :: use_sj_sedimentation_solver = .FALSE.
logical            :: dust_debug = .false.
logical            :: do_esm_dust_flux = .false. !If set to .true. a dust flux will be prepared for exchange with Ocean
integer            :: logunit

namelist /dust_nml/  dust_source_filename, dust_source_name, uthresh, coef_emis, use_sj_sedimentation_solver, &
                     dust_debug, do_esm_dust_flux

!---- version number -----
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------

contains

!#######################################################################
! this subroutine calculates tendencies for all dust tracers, and reports
! total fields, like total dust emission and settling
subroutine atmos_dust_sourcesink ( lon, lat, frac_land, pwt, dt, &
       zhalf, pfull, w10m, t, rh, tracer, dsinku, rdt, Time, is,ie,js,je, kbot)

  real, intent(in) :: lon(:,:), lat(:,:) ! geographical coordinates, units?
  real, intent(in) :: frac_land(:,:) ! fraction of land in the grid cell
  real, intent(in) :: w10m(:,:) ! 10-m wind speed, m/s
  real, intent(in) :: pwt(:,:,:) ! mass of each atmos layer, kg/m2
  real, intent(in) :: zhalf(:,:,:) ! z of half-layers, m(?)
  real, intent(in) :: pfull(:,:,:) ! pressure on layers, Pa
  real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degK
  real, intent(in) :: rh(:,:,:) ! relative humidity 
  real, intent(in) :: tracer(:,:,:,:) ! tracer concentrations
  real, intent(in) :: dsinku(:,:,:) ! dry deposition flux at the surface, for diag only
  real, intent(in) :: dt ! model timestep
  real, intent(inout) :: rdt(:,:,:,:) ! tendency of tracers, to be updated for dust tracers
  type(time_type), intent(in) :: Time ! current model time
  integer, intent(in) :: is, ie, js, je ! boundaries of physical window
  integer, intent(in), optional :: kbot(:,:) ! index of bottom level
  ! NOTE that operations are done on physics window; in particular the sizes
  ! of the arrays correspond to the sizes of physics window.

! ---- local vars
  real, dimension(size(tracer,1),size(tracer,2)) :: &
     source, &        ! source fraction
     all_dust_setl, & ! total dust sedimentation flux at the bottom of the atmos
     dust_emis, &     ! dust emission flux at the bottom of the atmos
     all_dust_emis    ! total dust emission flux at the bottom of the atmos
  real, dimension(size(tracer,1),size(tracer,2),size(tracer,3)) :: &
     dust_dt           ! calculated dust tendency

  integer :: i
  integer :: kd    ! vertical size of our arrays
  integer :: ndust ! atmos tracer number that corresponds to the current dust tracer
  logical :: used
  
  kd = size(tracer,3)
  
  ! initialize accumulated deposition and emission fields
  all_dust_emis(:,:) = 0.0
  all_dust_setl(:,:) = 0.0

  !----------- dust sources on local grid
  source(:,:)=0.0
  if (do_emission) then
     call interpolator(dust_source_interp, Time, source, trim(dust_source_name(1)), is, js)
  endif
  
  ! Send the dust source data to the diag_manager for output.
  if (id_dust_source > 0 ) &
       used = send_data ( id_dust_source, source, Time, is_in=is, js_in=js )

  do i = 1,n_dust_tracers
     ndust = dust_tracers(i)%tr
     ! calculate sources and sinks for each individual dust tracer
     call atmos_dust_sourcesink1(frac_land, pwt, dt, &
        dust_tracers(i)%dustden, dust_tracers(i)%dustref, dust_tracers(i)%frac_s, source, &
        pfull, w10m, t, rh, &
        tracer(:,:,:,ndust), dust_dt, dust_emis, &
        dust_tracers(i)%dust_setl(is:ie,js:je), dust_tracers(i)%dsetl_dtr(is:ie,js:je), &
        dust_tracers(i)%do_surf_exch, &
        is,ie,js,je, kbot)
     ! update dust tendencies
     rdt(:,:,:,ndust)=rdt(:,:,:,ndust)+dust_dt(:,:,:)
     
     ! Send the emission data to the diag_manager for output.
     if (dust_tracers(i)%id_dust_emis > 0 ) then
       used = send_data ( dust_tracers(i)%id_dust_emis, dust_emis, Time, is_in=is,js_in=js )
     endif
     ! Send the settling data to the diag_manager for output.
     if (dust_tracers(i)%id_dust_setl > 0 ) then
       used = send_data ( dust_tracers(i)%id_dust_setl, dust_tracers(i)%dust_setl(is:ie,js:je), Time, is_in=is,js_in=js )
     endif

     ! Accumulate total emission and deposition fluxes for output
     ! accumulate total dust deposition flux
     all_dust_setl(:,:) = all_dust_setl(:,:) &
          + dust_tracers(i)%dust_setl(is:ie,js:je) + pwt(:,:,kd)*dsinku(:,:,ndust) ! shouldn't kd be kbot?

!     all_dust_setl(:,:) = 1.0-frac_land(:,:)  !The exchanged flux of this becomes >1 at some points within ocean near shore!!
     
     if (id_dust_emis > 0 .or. id_emidust > 0) then
        ! accumulate total dust emission flux
        all_dust_emis(:,:) = all_dust_emis(:,:) + dust_emis(:,:) 
     endif
  enddo

  call atmos_dust_drydep_flux_set(all_dust_setl, is,ie,js,je)
 
  if (id_dust_ddep > 0) then
     used = send_data (id_dust_ddep, all_dust_setl(:,:), Time, is_in=is, js_in=js)
  endif
  if (id_dust_emis > 0) then
     used = send_data (id_dust_emis, all_dust_emis(:,:), Time, is_in=is, js_in=js)
  endif

  ! cmip variables
  if (id_drydust > 0) then
     used = send_data (id_drydust, all_dust_setl(:,:), Time, is_in=is, js_in=js)
  endif
  if (id_emidust > 0) then
     used = send_data (id_emidust, all_dust_emis(:,:), Time, is_in=is, js_in=js)
  endif

end subroutine atmos_dust_sourcesink


!#######################################################################
! Given dust properties and atmospheric variables, calculate dust tendencies
subroutine atmos_dust_sourcesink1 ( &
       frac_land, pwt, dt, &
       dustden, dustref, frac_s, source, &
       pfull, w10m, t, rh, &
       dust, dust_dt, dust_emis, dust_setl, dsetl_dtr, do_surf_exch, is,ie,js,je,kbot)

  real, intent(in),  dimension(:,:)   :: frac_land
  real, intent(in) :: dt ! model timestep
  real, intent(in) :: dustref ! effective radius of the dry dust particles, m
  real, intent(in) :: dustden ! density of dry dust particles, kg/m3
  real, intent(in) :: frac_s  ! fraction of source
  real, intent(in) :: source(:,:) ! dust source at the surface
  real, intent(in),  dimension(:,:)   :: w10m
  real, intent(out) :: dust_emis(:,:) ! dust emission
  real, intent(out) :: dust_setl(:,:) ! grav. sedimentation flux at the atmos bottom 
  real, intent(out) :: dsetl_dtr(:,:) ! derivative of dust_setl w.r.t. dust concentration
  real, intent(in),  dimension(:,:,:) :: pwt, dust
  real, intent(in),  dimension(:,:,:) :: pfull, t, rh
  real, intent(out), dimension(:,:,:) :: dust_dt
  integer, intent(in),  dimension(:,:), optional :: kbot
  integer, intent(in)  :: is, ie, js, je
  logical, intent(in)  :: do_surf_exch

  ! ---- local vars
  integer :: outunit, unit, ierr, io
  integer  i, j, k, id, jd, kd, kb
  real, dimension(size(dust,3)) :: vdep, setl, dz, air_dens, qn, qn1

  real, parameter :: mtcm = 100.  ! meter to cm
  real, parameter :: mtv  = 1.    ! factor conversion for mixing ratio of dust
  real, parameter :: ptmb = 0.01  ! pascal to mb

  real :: rhb, rcm
  real :: ratio_r, rho_wet_dust
  real :: rho_air
  real :: rwet

  id=size(dust,1); jd=size(dust,2); kd=size(dust,3)

  
  dust_emis(:,:) = 0.0
  dust_setl(:,:) = 0.0
  dsetl_dtr(:,:) = 0.0
  dust_dt(:,:,:) = 0.0

!----------- compute dust emission ------------
  where ( frac_land > 0.1 .and. w10m > u_ts ) &
      dust_emis = CH * frac_s * source * frac_land * w10m**2 * (w10m - u_ts)

  dust_dt(:,:,kd)=dust_dt(:,:,kd)+dust_emis(:,:)/pwt(:,:,kd)*mtv


  rcm=dustref*mtcm            ! Particles radius in centimeters
!------------------------------------------
!       Solve at the model TOP (layer plev-10)
!------------------------------------------
  do j=1,jd
  do i=1,id
     setl(:)=0.
     air_dens(:)=pfull(i,j,:)/t(i,j,:)/RDGAS
     dz(:) = pwt(i,j,:)/air_dens(:)
     if (present(kbot)) then
        kb=kbot(i,j)
     else
        kb=kd
     endif
     do k=1,kb
        rhb=amin1(0.99,rh(i,j,k))
        rhb=amax1(0.01,rhb)
!----------------------------------------------------------
!     Aerosol growth with relative humidity
!----------------------------------------------------------
        rwet         = dustref                     ! Add any particle growth here
        ratio_r      = (dustref/rwet)**3           ! Ratio dry over wet radius cubic power
        rho_wet_dust = ratio_r*dustden+(1.-ratio_r)*DENS_H2O ! Density of wet aerosol [kg/m3]
        vdep(k)      = sedimentation_velocity(t(i,j,k),pfull(i,j,k),rwet,rho_wet_dust) ! Settling velocity [m/s]
     enddo
     if (use_sj_sedimentation_solver) then
       qn(:)=dust(i,j,:)
       qn1(1)=qn(1)*dz(1)/(dz(1)+dt*vdep(1))
       do k=2,kb 
         qn1(k)=(qn(k)*dz(k)+dt*qn1(k-1)*vdep(k-1)*air_dens(k-1)/air_dens(k))/(dz(k)+dt*vdep(k))
       enddo
       dust_dt(i,j,:)=dust_dt(i,j,:)+(qn1(:)-qn(:))/dt
       setl(kb) = qn1(kb)*air_dens(kb)/mtv*vdep(kb) !! ter, 2016-05-13 vdep(k) is outside of loop, so k > kb.  Switched k to kb 

!---> h1g, 2016-04-05
       if( dust_debug ) then
!$OMP CRITICAL

     !  if (mpp_pe()==mpp_root_pe()) then
     !      write(logunit,'("DUST ",2i5," qn(kb)=",e12.4," qn1(kb)=",e12.4," vdep(kb)=",e12.4," air_dens(kb)=",e12.4," dz(kb)=",e12.4," dt=",e12.4," dust_dt=",e12.4," setl=",e12.4)')  &
     !               i,j,qn(kb),qn1(kb),vdep(kb),air_dens(kb),dz(kb),dt,dust_dt(i,j,kb),setl(kb)
     !  endif
!$OMP END CRITICAL  
       end if ! dust_debug
!<--- h1g, 2016-04-05

     else
       do k=1,kb
        if (dust(i,j,k) > 0.0) then
          setl(k)=dust(i,j,k)*air_dens(k)/mtv*vdep(k)    ! settling flux [kg/m2/s]
        endif
       enddo
       dust_dt(i,j,1)=dust_dt(i,j,1)-setl(1)/pwt(i,j,1)*mtv
       dust_dt(i,j,2:kb)=dust_dt(i,j,2:kb) &
          + ( setl(1:kb-1) - setl(2:kb) )/pwt(i,j,2:kb)*mtv
     endif
     dust_setl(i,j)  = setl(kb)          ! at the bottom of the atmos
     rho_air = air_dens(kb)           ! rho_air is not defined
     dsetl_dtr(i,j) = rho_air/mtv*vdep(kb) ! derivative of settling flux w.r.t tracer conc
     if (do_surf_exch) &
       setl(kb) = setl(kb)*(1-frac_land(i,j)) ! settlement tendency in the  
       ! near-surface layer over land will be handled by flux exchange
  enddo
  enddo 
end subroutine atmos_dust_sourcesink1


!#######################################################################
! calculates the vertical velocity of dust settling
elemental real function sedimentation_velocity(T,p,rwet,rho_wet_dust) result(vdep)
   real, intent(in) :: T            ! air temperature, deg K
   real, intent(in) :: p            ! pressure, Pa
   real, intent(in) :: rwet         ! radius of dust particles, m
   real, intent(in) :: rho_wet_dust ! density of dust particles, kg/m3
 
   real :: viscosity, free_path, C_c
   viscosity = 1.458E-6 * T**1.5/(T+110.4)     ! Dynamic viscosity
   free_path = 6.6e-8*T/293.15*(PSTD_MKS/p)
   C_c = 1.0 + free_path/rwet * &              ! Slip correction [none]
               (1.257+0.4*exp(-1.1*rwet/free_path))
   vdep = 2./9.*C_c*GRAV*rho_wet_dust*rwet**2/viscosity  ! Settling velocity [m/s]
end function sedimentation_velocity


!######################################################################
! given a tracer index, returns TRUE if this is one of dust tracers
function is_dust_tracer(tr) result(ret)
  logical :: ret
  integer, intent(in) :: tr

  integer :: i
  ret = .FALSE.
  do i = 1,n_dust_tracers
     ret = (dust_tracers(i)%tr==tr)
     if (ret) return
  enddo
end function 


!######################################################################
! given a tracer index, returns TRUE if this is one of dust tracers
function dust_has_surf_setl_flux(tr) result(ret)
  logical :: ret
  integer, intent(in) :: tr

  integer :: i
  ret = .FALSE.
  do i = 1,n_dust_tracers
     if (dust_tracers(i)%tr==tr) then
        ret=dust_tracers(i)%do_surf_exch
        return
     endif
  enddo   
end function 


!######################################################################
subroutine get_dust_surf_setl_flux(tr, dust_setl, dsetl_dtr)
  integer, intent(in)  :: tr ! tracer index
  real,    intent(out) :: dust_setl(:,:) ! sedimentation flux, on compute domain
  real,    intent(out) :: dsetl_dtr(:,:) ! derivative of dust_setl w.r.t. dust concentration
  
  integer :: i
  dust_setl(:,:) = 0.0
  dsetl_dtr(:,:) = 0.0
  do i = 1,n_dust_tracers
     if (dust_tracers(i)%tr/=tr) cycle
     dust_setl(:,:) = dust_tracers(i)%dust_setl(:,:)
     dsetl_dtr(:,:) = dust_tracers(i)%dsetl_dtr(:,:)
     return
  enddo
end subroutine 


!#######################################################################

!<SUBROUTINE NAME="atmos_dust_init">
!<OVERVIEW>
! The constructor routine for the dust module.
!</OVERVIEW>
subroutine atmos_dust_init (lonb, latb, axes, Time, mask)
  real,             intent(in) :: lonb(:,:), latb(:,:) ! grid cell boundaries
  type(time_type),  intent(in) :: Time                 ! model time
  integer,          intent(in) :: axes(4)              ! diagnostic axes
  real, optional,   intent(in) :: mask(:,:,:)

  ! ---- local vars
  integer :: outunit, unit, ierr, io
  integer :: n_atm_tracers ! number of prognostic atmos tracers
  integer :: tr ! atmos tracer iterator
  integer :: i  ! running index of dust tracers
  character(32)  :: tr_name    ! tracer name
  character(256) :: longname   ! long and meaningful tracer name
  character(32)  :: method     ! method string for parameter retrieval (not used)
  character(256) :: parameters ! parameter string for dust tracer
  real    :: value ! temporary storage for parsing input
  real    :: frac_s_sum ! sum of source fractions frac_s
  
  if (module_is_initialized) return

  call write_version_number (version, tagname)
  logunit = stdlog()
  outunit = stdout()

  ! read namelist.
  if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=dust_nml, iostat=io)
    ierr = check_nml_error(io,'dust_nml')
#else
    unit =  open_namelist_file ( )
    ierr=1; do while (ierr /= 0)
       read (unit, nml=dust_nml, iostat=io, end=10)
       ierr = check_nml_error(io, 'dust_nml')
    end do
10  call close_file (unit)
#endif
  endif
 
  ! write namelist to the log file
  if (mpp_pe() == mpp_root_pe()) then
     write (logunit, nml=dust_nml)
  endif

  if (uthresh .le. -990) then
    u_ts = 0.
  else
    u_ts=uthresh
  endif
  if (coef_emis .le. -990) then
    ch = 1.0e-10
  else
    ch = coef_emis
  endif

  ! find out if there are any dust tracers
  call get_number_tracers(MODEL_ATMOS, num_prog=n_atm_tracers)
  n_dust_tracers = 0
  do tr = 1, n_atm_tracers
     call get_tracer_names(MODEL_ATMOS,tr,tr_name)
     if (lowercase(tr_name(1:4))=='dust') then
        n_dust_tracers = n_dust_tracers + 1
     endif
  enddo
  
  ! log the number of dust tracers
  call error_mesg('atmos_dust_init','number of atmos dust tracers ='//string(n_dust_tracers), &
                  NOTE)

  ! check if any dust is present in the model and set up the flag
  if (n_dust_tracers == 0) return ! nothing to do

  ! get the size of our compute domain, for settl flux storage
  
  ! fill the tracer parameters
  allocate(dust_tracers(n_dust_tracers))
  i = 0
  ierr = 0 
  do tr = 1, n_atm_tracers
     call get_tracer_names(MODEL_ATMOS,tr,name=tr_name,longname=longname)
     if (lowercase(tr_name(1:4)).ne.'dust') cycle ! this is not dust, we are not interested  

     i = i+1
     dust_tracers(i)%name = tr_name
     dust_tracers(i)%tr   = tr
     call set_tracer_atts(MODEL_ATMOS,tr_name,longname,'mmr')
     ! allocate space to store dust sedimentation flux for exchange with land.
     ! sizes of lonb and latb are used to get the size of the compute domain  
     allocate(dust_tracers(i)%dust_setl(size(lonb)-1,size(latb)-1))
     allocate(dust_tracers(i)%dsetl_dtr(size(lonb)-1,size(latb)-1))

     method = ''; parameters = ''
     if(query_method('parameters', MODEL_ATMOS, tr, method, parameters)) then
        call parse_and_check(parameters, tr_name, 'ra',      dust_tracers(i)%ra,      ierr)
        call parse_and_check(parameters, tr_name, 'rb',      dust_tracers(i)%rb,      ierr)
        call parse_and_check(parameters, tr_name, 'dustref', dust_tracers(i)%dustref, ierr)
        call parse_and_check(parameters, tr_name, 'dustden', dust_tracers(i)%dustden, ierr)
     else
        call error_mesg('atmos_dust_init',&
          '"parameters" line is missing from the field table for dust tracer "'//trim(tr_name)//'"', NOTE)
        ierr = ierr+1
     endif
     
     if(query_method('emission', MODEL_ATMOS, tr, method, parameters)) then
        do_emission = (trim(method) == 'prescribed')
        if (do_emission) & 
           call parse_and_check(parameters, tr_name, 'source_fraction',  dust_tracers(i)%frac_s, ierr)
     else
        call error_mesg('atmos_dust_init',&
          '"emission" line is missing from the field table for dust tracer "'//trim(tr_name)//'"', NOTE)
        ierr = ierr+1
     endif

     if (query_method('dry_deposition', MODEL_ATMOS, tr, method, parameters)) then
        dust_tracers(i)%do_surf_exch = (index(method,'lm3')>0)
        ! we do exchange with surface if "lm3" appears in the dry deposition
        ! method in the field table
        if (dust_tracers(i)%do_surf_exch) then
           ! check that land has the tracer for the exchange with the atmos
           if ( get_tracer_index ( MODEL_LAND, tr_name ) == NO_TRACER ) &
              call error_mesg ('atmos_dist_init',&
                 'dry deposition method "lm3" requires presence of tracer "'&
                 //trim(tr_name)//'" in both atmos and land tracer tables.', FATAL)
        endif
     endif

     ! Register a diagnostic field : total emission of dust
     dust_tracers(i)%id_dust_emis = register_diag_field ( module_name,     &
                     trim(dust_tracers(i)%name)//'_emis', axes(1:2),Time,  &
                     trim(dust_tracers(i)%name)//'_emis', 'kg/m2/s',       &
                     missing_value=-999.  )
     ! Register a diagnostic field : total settling of dust
     dust_tracers(i)%id_dust_setl = register_diag_field ( module_name,     &
                     trim(dust_tracers(i)%name)//'_setl', axes(1:2),Time,  &
                     trim(dust_tracers(i)%name)//'_setl', 'kg/m2/s',       &
                     missing_value=-999.  )
  enddo  
  ! print out information about dust tracers
  if (mpp_pe()==mpp_root_pe()) then
     call print_table(logunit)
     call print_table(outunit)
  endif
  ! stop the models is there were any errors reading dust parameters
  if (ierr>0) &
     call error_mesg('atmos_dust_init', &
         'ERRORS detected reading one or more parameters of dust tracers; '// &
         'look for NOTES from atmos_dust_init above',&
         FATAL)
!  normalize source fraction
!  NOTE: normalization is commented out because apparently in the original code 
!        frac_s coefficients were not normalized:
!        0.05+0.1125+0.225+0.225+0.225 = 0.8375
!  frac_s_sum = sum(dust_tracers(1:n_dust_tracers)%frac_s)
!  if (frac_s_sum <= 0) &
!     call error_mesg('atmos_dust_init', &
!         'sum of per-bin source fractions frac_s must be positive; it is '//string(frac_s_sum), &
!         FATAL)
!  dust_tracers(1:n_dust_tracers)%frac_s = dust_tracers(1:n_dust_tracers)%frac_s/frac_s_sum
  
  id_dust_source  = register_diag_field ( module_name,             &
                   'DU_source', axes(1:2), Time,                   &
                   'DU_source', 'none')
  id_dust_ddep = register_diag_field ( module_name, &
      'dust_ddep', axes(1:2), Time, &
      'total dry deposition and settling of dust', 'kg/m2/s')

  id_dust_emis = register_diag_field ( module_name, &
      'dust_emis', axes(1:2), Time, &
      'total emission of dust', 'kg/m2/s')

  ! register cmip variables
  id_drydust = register_cmip_diag_field_2d ( module_name, 'drydust', Time, &
                              'Dry Deposition Rate of Dust', 'kg m-2 s-1', &
                standard_name='tendency_of_atmosphere_mass_content_of_dust_dry_aerosol_particles_due_to_dry_deposition')

  id_emidust = register_cmip_diag_field_2d ( module_name, 'emidust', Time, &
                              'Total Emission Rate of Dust', 'kg m-2 s-1', &
                standard_name='tendency_of_atmosphere_mass_content_of_dust_dry_aerosol_particles_due_to_emission')

 
  if (do_emission) then
     call interpolator_init (dust_source_interp, trim(dust_source_filename),  &
                          lonb, latb, &
                          data_out_of_bounds=  (/CONSTANT/), &
                          data_names = dust_source_name, &
                          vert_interp=(/INTERP_WEIGHTED_P/) )
  endif

  !Allocate the array to contain the total dust flux for ESM
  allocate(dry_dep_lith_dust_flux(size(lonb)-1,size(latb)-1)); dry_dep_lith_dust_flux=0.0
  allocate(wet_dep_lith_dust_flux(size(lonb)-1,size(latb)-1)); wet_dep_lith_dust_flux=0.0


  do_dust = .TRUE.
  module_is_initialized = .TRUE.

 end subroutine atmos_dust_init
!</SUBROUTINE>
!#######################################################################

!<SUBROUTINE NAME ="atmos_dust_flux_init">

!<OVERVIEW>
! Subroutine to initialize the dust flux
!</OVERVIEW>

 subroutine atmos_dust_flux_init

   use tracer_manager_mod,    only : get_tracer_index, tracer_manager_init
   use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux

   !
   !-----------------------------------------------------------------------
   !     arguments
   !-----------------------------------------------------------------------
   !

   !
   !-----------------------------------------------------------------------
   !     local variables
   !-----------------------------------------------------------------------
   !

   integer :: ind

   !
   !-----------------------------------------------------------------------
   !     local parameters
   !-----------------------------------------------------------------------
   !

   character(len=64), parameter    :: mod_name = 'atmos_dust'
   character(len=64), parameter    :: sub_name = 'atmos_dust_flux_init'
   character(len=256), parameter   :: error_header =                               &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: warn_header =                                &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

   integer :: outunit
   outunit = stdout()
   
   if(do_esm_dust_flux) then
   ! find out if there is a dust tracer (to be used to pass dust fluxes to Ocean BGC)
   call tracer_manager_init      
   ind = get_tracer_index(MODEL_ATMOS,'dust1')
   if (ind > 0) then
      write (outunit,*) trim(note_header), ' dust1 was initialized as tracer number ', ind
   else
      write (outunit,*) trim(note_header), ' dust1 was not found ', ind
      return
   endif

   !
   !       Initialize coupler flux for the total dust flux to ocean.
   !       Use "dust1" tracer if it exists to avoid creating a dummy tracer for this purpose.
   !       This does not affect the existing "dust1" tracer unless "dust1" should need to have
   !       its own coupler fluxes to be exchanged with OCN in the future. In that case we need
   !       to define an additional prog.tracer in ATM to associate the ESM-wanted "lith" fluxes with.
   !

   if (ind > 0) then
      ind_dry_dep_lith_dust_flux = aof_set_coupler_flux('dry_dep_lith', &
           flux_type = 'air_sea_deposition', implementation = 'dry',    &
           atm_tr_index = ind,                                          &
           mol_wt = 1.0, param = (/ 1.0,1.0 /),                         &
           caller = trim(mod_name) // '(' // trim(sub_name) // ')')

      ind_wet_dep_lith_dust_flux = aof_set_coupler_flux('wet_dep_lith', &
           flux_type = 'air_sea_deposition', implementation = 'wet',    &
           atm_tr_index = ind,                                          &
           mol_wt = 1.0, param = (/ 1.0,1.0 /),                         &
           caller = trim(mod_name) // '(' // trim(sub_name) // ')')
   endif
   endif !if(do_esm_dust_flux)

 end subroutine atmos_dust_flux_init
!</SUBROUTINE>
!#######################################################################

!<SUBROUTINE NAME ="atmos_dust_gather_data">
!<OVERVIEW>
!  A subroutine to gather fields needed for calculating the dust flux for ESM
!</OVERVIEW>
!

subroutine atmos_dust_gather_data (gas_fields, tr_bot)

use coupler_types_mod, only: coupler_2d_bc_type, ind_pcair

type(coupler_2d_bc_type), intent(inout) :: gas_fields
real, dimension(:,:,:), intent(in)      :: tr_bot


if (ind_dry_dep_lith_dust_flux .gt. 0) then
  gas_fields%bc(ind_dry_dep_lith_dust_flux)%field(ind_pcair)%values(:,:) = -dry_dep_lith_dust_flux(:,:)!sign flip
endif

if (ind_wet_dep_lith_dust_flux .gt. 0) then
  gas_fields%bc(ind_wet_dep_lith_dust_flux)%field(ind_pcair)%values(:,:) = wet_dep_lith_dust_flux(:,:)
endif

end subroutine atmos_dust_gather_data
!</SUBROUTINE >

!#######################################################################

!#######################################################################
! read the specified parameter from the input string, and report a
! non-fatal error if something goes wrong
subroutine parse_and_check(parameters, tr_name, par_name, value, nerrors)
  character(*), intent(in) :: &
     parameters, &              ! string to parse
     tr_name,    &              ! tracer name, for reporting only
     par_name                   ! name of the parameter
  real, intent(inout)  :: value ! result of the parsing
  integer, intent(inout) :: nerrors ! error flag

  real :: tmp ! buffer for input value
  
  if ( parse(parameters, par_name, tmp) > 0 ) then 
     value = tmp
  else
     nerrors = nerrors + 1
     call error_mesg('atmos_dust_init', &
       'ERROR reading parameter "'//trim(par_name)//'" for tracer "'//trim(tr_name)//'"',&
       NOTE)
  endif
end subroutine


!######################################################################
subroutine print_table(unit)
   integer, intent(in) :: unit
   
   integer :: i

   write(unit,'(x,121("-"))')
   write(unit,'(3x,99(x,a16))')'dust tr. name','atm. tr. number','ra','rb', &
       'dustref','dustden','frac_s','do_surf_exch'
   write(unit,'(x,121("-"))')
   do i = 1,n_dust_tracers
      write(unit,'(x,i2,x,a16,99(x,g16.6))')&
         i, trim(dust_tracers(i)%name), dust_tracers(i)%tr, &
         dust_tracers(i)%ra, dust_tracers(i)%rb, &
         dust_tracers(i)%dustref, dust_tracers(i)%dustden, &
         dust_tracers(i)%frac_s, dust_tracers(i)%do_surf_exch
   enddo
   write(unit,'(x,121("-"))')   
end subroutine 


!######################################################################
subroutine atmos_dust_time_vary (Time)
  type(time_type), intent(in) :: Time

  if (do_emission) &
     call obtain_interpolator_time_slices (dust_source_interp, Time)

end subroutine atmos_dust_time_vary 


!######################################################################
subroutine atmos_dust_endts              
  if (do_emission) &
      call unset_interpolator_time_flag (dust_source_interp)
end subroutine atmos_dust_endts 

!######################################################################

subroutine atmos_dust_wetdep_flux_set(array, is,ie,js,je)
  real, dimension(is:ie,js:je), intent(in) :: array
  integer,              intent(in) :: is,ie,js,je
  if (n_dust_tracers == 0) return ! nothing to do
  wet_dep_lith_dust_flux(is:ie,js:je) = array(is:ie,js:je)
end subroutine atmos_dust_wetdep_flux_set

subroutine atmos_dust_drydep_flux_set(array, is,ie,js,je)
  real, dimension(is:ie,js:je), intent(in) :: array
  integer,              intent(in) :: is,ie,js,je
  if (n_dust_tracers == 0) return ! nothing to do
  dry_dep_lith_dust_flux(is:ie,js:je) = array(is:ie,js:je)
end subroutine atmos_dust_drydep_flux_set


!#######################################################################
!<SUBROUTINE NAME="atmos_dust_end">
!<OVERVIEW>
!  The destructor routine for the dust module.
!</OVERVIEW>
 subroutine atmos_dust_end
 
    integer :: i
    if (do_emission) &
         call interpolator_end (dust_source_interp)
    module_is_initialized = .FALSE.
    if(.not. do_dust) return
    do_dust = .FALSE.
    do i = 1,n_dust_tracers
       deallocate(dust_tracers(i)%dust_setl)
       deallocate(dust_tracers(i)%dsetl_dtr)
    enddo
    deallocate(dust_tracers)
    deallocate(dry_dep_lith_dust_flux)
    deallocate(wet_dep_lith_dust_flux)
 end subroutine atmos_dust_end
!</SUBROUTINE>

end module atmos_dust_mod
