module atmos_soa_mod
! <DESCRIPTION>
!   This module is an implementation of Secondary organic aerosols (SOA)
!   from anthropogenic activities, and is based on Tie et al. (JGR, 2003).
!   The only souce of SOA is due to the oxydation of C4H10 by OH.
!   The concentrations of these 2 gas species are read as input.
! </DESCRIPTION>
! <WARNING>
!  To save space only the actual month of input files are kept in memory. 
!  This implies that the "atmos_SOA_init" should be executed at the begining 
!  of each month. In other words, the script should not run more than 1 month
!  without a restart.
! </WARNING>
! <CONTACT EMAIL="Paul.Ginoux@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use mpp_mod, only: input_nml_file 
use                    fms_mod, only : file_exist,              &
                                       write_version_number,    &
                                       mpp_pe,                  &
                                       mpp_root_pE,             &
                                       close_file,              &
                                       stdlog,                  &
                                       check_nml_error, error_mesg, &
                                       open_namelist_file, FATAL, NOTE
use           time_manager_mod, only : time_type,               &
                                       days_in_month, &
                                       set_date, set_time, print_date, get_date, &
                                       operator(>), operator(+), operator(-)
use           diag_manager_mod, only : send_data,               &
                                       register_diag_field,     &
                                       register_static_field,   &
                                       get_base_time
use        atmos_cmip_diag_mod, only : register_cmip_diag_field_2d
use         tracer_manager_mod, only : get_tracer_index,        &
                                       set_tracer_atts
use          field_manager_mod, only : MODEL_ATMOS
use              constants_mod, only : PI, GRAV, RDGAS, WTMAIR, AVOGNO
use           interpolator_mod, only:  interpolate_type,  &
                                       interpolator_init, &
                                       obtain_interpolator_time_slices,&
                                       unset_interpolator_time_flag, &
                                       interpolator, interpolator_end, &
                                       CONSTANT, INTERP_WEIGHTED_P

implicit none

private
!-----------------------------------------------------------------------
!----- interfaces -------
!
public  atmos_SOA_init, atmos_SOA_end, atmos_SOA_chem, &
        atmos_SOA_time_vary, atmos_soa_endts

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

integer :: nSOA = 0  ! tracer number for Secondary Organic Aerosol 

!--- identification numbers for  diagnostic fields and axes ----
integer ::   id_OH_conc            = 0
integer ::   id_C4H10_conc         = 0
integer ::   id_SOA_chem           = 0
integer ::   id_SOA_chem_col       = 0
integer ::   id_chepsoa            = 0 ! cmip field (same as SOA_chem_col)
integer ::   id_SOA_isoprene       = 0
integer ::   id_SOA_terpene        = 0
integer ::   id_SOA_biogenic       = 0

type(interpolate_type),save         ::  gas_conc_interp
type(interpolate_type),save         ::  isoprene_interp
type(interpolate_type),save         ::  terpene_interp

type(time_type) :: isoprene_offset
type(time_type) :: terpene_offset
type(time_type) :: isoprene_entry
type(time_type) :: terpene_entry
logical :: isoprene_negative_offset
logical :: terpene_negative_offset
integer :: isoprene_time_series_type
integer :: terpene_time_series_type
type(time_type) :: isoprene_time
type(time_type) :: terpene_time

! Initial calendar time for model
type(time_type) :: model_init_time

real, parameter       :: wtm_C = 12.
real, parameter       :: wtm_C4H10 = 58.
integer, parameter    :: carbon_per_isoprene = 5
integer, parameter    :: carbon_per_terpene = 10
real, parameter       :: cm2_per_m2 = 1.e4
real, parameter       :: kg_per_g = 1.e-3
real, parameter       :: om_oc_ratio = 1.5
! Factors to convert from molecules(BVOC)/cm2/s to kg(OM)/m2/s
real, parameter       :: isoprene_factor = cm2_per_m2 / AVOGNO * wtm_C * kg_per_g * carbon_per_isoprene * om_oc_ratio
real, parameter       :: terpene_factor = cm2_per_m2 / AVOGNO * wtm_C * kg_per_g * carbon_per_terpene * om_oc_ratio
integer, parameter    :: TS_CONSTANT=1, TS_FIXED=2, TS_VARYING=3

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------

character(len=32)  :: gas_conc_filename = 'gas_conc_3D.nc'
character(len=32)  :: isoprene_filename = 'biogenic_emis.nc'
character(len=32)  :: terpene_filename  = 'biogenic_emis.nc'
character(len=32), dimension(2) :: gas_conc_name
data gas_conc_name/'OH','C4H10'/
character(len=32), dimension(1) :: isoprene_input_name = (/ 'isoprene' /)
character(len=32), dimension(1) :: terpene_input_name = (/ 'terpenes' /)
integer, dimension(6) :: isoprene_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
integer, dimension(6) :: terpene_dataset_entry   = (/ 1, 1, 1, 0, 0, 0 /)
character(len=80)     :: isoprene_source = ' '
character(len=80)     :: terpene_source = ' '
character(len=80)     :: isoprene_time_dependency_type = 'constant'
character(len=80)     :: terpene_time_dependency_type  = 'constant'
real                  :: isoprene_SOA_yield = 0.05
real                  :: terpene_SOA_yield  = 0.05
namelist /secondary_organics_nml/ gas_conc_filename, gas_conc_name, &
                                  isoprene_source, &
                                  isoprene_filename, &
                                  isoprene_input_name, &
                                  isoprene_SOA_yield, &
                                  isoprene_time_dependency_type, &
                                  isoprene_dataset_entry, &
                                  terpene_source, &
                                  terpene_filename, &
                                  terpene_input_name, &
                                  terpene_SOA_yield, &
                                  terpene_time_dependency_type, &
                                  terpene_dataset_entry

logical :: module_is_initialized=.FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
!-----------------------------------------------------------------------

contains


!#######################################################################

!<SUBROUTINE NAME="atmos_SOA_init">
!<OVERVIEW>
! The constructor routine for the soa module.
!</OVERVIEW>
subroutine atmos_SOA_init ( lonb, latb, nlev, axes, Time, mask)
!-----------------------------------------------------------------------
real,             intent(in), dimension(:,:)        :: lonb, latb
integer,          intent(in)                        :: nlev
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
real, intent(in), dimension(:,:,:), optional        :: mask
character(len=7), parameter :: mod_name = 'tracers'
!
!-----------------------------------------------------------------------
!
      integer  unit,io,ierr, logunit
      character(len=3) :: SOA_tracer
!
      data SOA_tracer/'SOA'/

!
      if (module_is_initialized) return
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=secondary_organics_nml, iostat=io)
        ierr = check_nml_error(io,'secondary_organics_nml')
#else
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=secondary_organics_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'secondary_organics_nml')
        end do
10      call close_file (unit)
#endif
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit=stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=secondary_organics_nml)

!----- set initial value of soa ------------

      nSOA = get_tracer_index(MODEL_ATMOS,'SOA')
      if (nSOA > 0) then
         call set_tracer_atts(MODEL_ATMOS,'SOA','SOA','mmr')
         if (nSOA > 0 .and. mpp_pe() == mpp_root_pe()) &
                 write (*,30) SOA_tracer,nsoa
         if (nSOA > 0 .and. mpp_pe() == mpp_root_pe()) &
                 write (logunit,30) SOA_tracer,nsoa
      endif


  30   format (A,' was initialized as tracer number ',i2)

     call interpolator_init (gas_conc_interp, trim(gas_conc_filename),  &
                             lonb, latb,&        
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = gas_conc_name, & 
                             vert_interp=(/INTERP_WEIGHTED_P/) )



!----------------------------------------------------------------------
!    define the model base time  (defined in diag_table)
!----------------------------------------------------------------------
        model_init_time = get_base_time()

   if ( trim(isoprene_source) .ne. ' ') then

        isoprene_offset = set_time (0,0)
        isoprene_entry = set_time (0,0)
        isoprene_negative_offset = .false.
        isoprene_time_series_type = TS_CONSTANT

!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(isoprene_time_dependency_type) == 'constant' ) then
        isoprene_time_series_type = TS_CONSTANT
        isoprene_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'isoprene is constant in atmos_soa module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for isoprene is selected.
!---------------------------------------------------------------------
      else if (trim(isoprene_time_dependency_type) == 'time_varying') then
        isoprene_time_series_type = TS_VARYING
        if (isoprene_dataset_entry(1) == 1 .and. &
            isoprene_dataset_entry(2) == 1 .and. &
            isoprene_dataset_entry(3) == 1 .and. &
            isoprene_dataset_entry(4) == 0 .and. &
            isoprene_dataset_entry(5) == 0 .and. &
            isoprene_dataset_entry(6) == 0 ) then
          isoprene_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to isoprene_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          isoprene_entry  = set_date (isoprene_dataset_entry(1), &
                                  isoprene_dataset_entry(2), &
                                  isoprene_dataset_entry(3), &
                                  isoprene_dataset_entry(4), &
                                  isoprene_dataset_entry(5), &
                                  isoprene_dataset_entry(6))
        endif
        call print_date (isoprene_entry , str= &
          'Data from isoprene timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        isoprene_offset = isoprene_entry - model_init_time
        if (model_init_time > isoprene_entry) then
          isoprene_negative_offset = .true.
        else
          isoprene_negative_offset = .false.
        endif
      else if (trim(isoprene_time_dependency_type) == 'fixed_year') then
        isoprene_time_series_type = TS_FIXED
        if (isoprene_dataset_entry(1) == 1 .and. &
            isoprene_dataset_entry(2) == 1 .and. &
            isoprene_dataset_entry(3) == 1 .and. &
            isoprene_dataset_entry(4) == 0 .and. &
            isoprene_dataset_entry(5) == 0 .and. &
            isoprene_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_soa_mod', &
            'must set isoprene_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to isoprene_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        isoprene_entry  = set_date (isoprene_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_soa_mod', &
           'isoprene is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'isoprene correspond to year :', &
                    isoprene_dataset_entry(1)
        endif
     endif
     call interpolator_init (isoprene_interp,           &
       trim(isoprene_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
       data_names=isoprene_input_name,vert_interp=(/INTERP_WEIGHTED_P/))
   endif


   if ( trim(terpene_source) .ne. ' ') then

        terpene_offset = set_time (0,0)
        terpene_entry = set_time (0,0)
        terpene_negative_offset = .false.
        terpene_time_series_type = TS_CONSTANT

!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(terpene_time_dependency_type) == 'constant' ) then
        terpene_time_series_type = TS_CONSTANT
        terpene_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'terpene is constant in atmos_soa module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for terpene is selected.
!---------------------------------------------------------------------
      else if (trim(terpene_time_dependency_type) == 'time_varying') then
        terpene_time_series_type = TS_VARYING
        if (terpene_dataset_entry(1) == 1 .and. &
            terpene_dataset_entry(2) == 1 .and. &
            terpene_dataset_entry(3) == 1 .and. &
            terpene_dataset_entry(4) == 0 .and. &
            terpene_dataset_entry(5) == 0 .and. &
            terpene_dataset_entry(6) == 0 ) then
          terpene_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to terpene_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          terpene_entry  = set_date (terpene_dataset_entry(1), &
                                  terpene_dataset_entry(2), &
                                  terpene_dataset_entry(3), &
                                  terpene_dataset_entry(4), &
                                  terpene_dataset_entry(5), &
                                  terpene_dataset_entry(6))
        endif
        call print_date (terpene_entry , str= &
          'Data from terpene timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        terpene_offset = terpene_entry - model_init_time
        if (model_init_time > terpene_entry) then
          terpene_negative_offset = .true.
        else
          terpene_negative_offset = .false.
        endif
      else if (trim(terpene_time_dependency_type) == 'fixed_year') then
        terpene_time_series_type = TS_FIXED
        if (terpene_dataset_entry(1) == 1 .and. &
            terpene_dataset_entry(2) == 1 .and. &
            terpene_dataset_entry(3) == 1 .and. &
            terpene_dataset_entry(4) == 0 .and. &
            terpene_dataset_entry(5) == 0 .and. &
            terpene_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_soa_mod', &
            'must set terpene_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to terpene_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        terpene_entry  = set_date (terpene_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_soa_mod', &
           'terpene is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'terpene correspond to year :', &
                    terpene_dataset_entry(1)
        endif
     endif
     call interpolator_init (terpene_interp,           &
       trim(terpene_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
       data_names=terpene_input_name,vert_interp=(/INTERP_WEIGHTED_P/))
   endif




      if (id_OH_conc .eq. 0 ) &
        id_OH_conc    = register_diag_field ( mod_name,           &
                      'OH_SOA_conc',axes(1:3),Time,                        &
                      'Hydroxyl radical concentration',           &
                      'molec.cm-3')

      id_C4H10_conc    = register_diag_field ( mod_name,           &
                      'C4H10_mmr',axes(1:3),Time,                        &
                      'nButane concentration',           &
                      'mmr')

      id_SOA_chem    = register_diag_field ( mod_name,       &
                      'SOA_chem',axes(1:3),Time,            &
                      'SOA production by C4H10 + OH',        &
                      'kg/m2/s')

      id_SOA_chem_col= register_diag_field ( mod_name,       &
                      'SOA_chem_col',axes(1:2),Time,            &
                      'column SOA production by C4H10 + OH',        &
                      'kg/m2/s')

      id_chepsoa = register_cmip_diag_field_2d ( mod_name, 'chepsoa', Time, &
                       'Chemical Production of Dry Aerosol Secondary Organic Matter', 'kg m-2 s-1', &
                      standard_name='tendency_of_atmosphere_mass_content_of_secondary_particulate_organic_matter_dry_aerosol_particles_due_to_net_chemical_production')
                      
      id_SOA_isoprene  = register_diag_field ( mod_name,       &
                        'SOA_isoprene',axes(1:2),Time,            &
                        'SOA pseudo-emission from isoprene',        &
                        'kg/m2/s')

      id_SOA_terpene   = register_diag_field ( mod_name,       &
                        'SOA_terpene',axes(1:2),Time,            &
                        'SOA pseudo-emission from terpene',        &
                        'kg/m2/s')

      id_SOA_biogenic  = register_diag_field ( mod_name,       &
                        'SOA_biogenic',axes(1:2),Time,            &
                        'total SOA pseudo-emission from biogenic VOCs',        &
                        'kg/m2/s')

      call write_version_number (version, tagname)

      module_is_initialized = .TRUE.

!-----------------------------------------------------------------------
end subroutine atmos_SOA_init




!#####################################################################

subroutine atmos_SOA_time_vary (Time)

type(time_type), intent(in) :: Time

    integer ::  yr, dum, mo_yr, mo, dy, hr, mn, sc, dayspmn


      call obtain_interpolator_time_slices (gas_conc_interp, Time)

   if ( trim(isoprene_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the isoprene data set from which data is to be 
!    taken. if isoprene is not time-varying, it is simply Time.
!---------------------------------------------------------------------
     if(isoprene_time_series_type .eq. TS_VARYING) then
       if (isoprene_negative_offset) then
         isoprene_time = Time - isoprene_offset
       else
         isoprene_time = Time + isoprene_offset
       endif
     else 
       if(isoprene_time_series_type .eq. TS_FIXED ) then
         call get_date (isoprene_entry, yr, dum,dum,dum,dum,dum)
         call get_date (Time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(isoprene_entry)
           if (dayspmn /= 29) then
             isoprene_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             isoprene_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           isoprene_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         isoprene_time = Time
       endif
     endif
     call obtain_interpolator_time_slices   &
                                     (isoprene_interp, isoprene_time)
   endif

   if ( trim(terpene_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the terpene data set from which data is to be 
!    taken. if terpene is not time-varying, it is simply Time.
!---------------------------------------------------------------------
     if(terpene_time_series_type .eq. TS_VARYING) then
       if (terpene_negative_offset) then
         terpene_time = Time - terpene_offset
       else
         terpene_time = Time + terpene_offset
       endif
     else 
       if(terpene_time_series_type .eq. TS_FIXED ) then
         call get_date (terpene_entry, yr, dum,dum,dum,dum,dum)
         call get_date (Time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(terpene_entry)
           if (dayspmn /= 29) then
             terpene_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             terpene_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           terpene_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         terpene_time = Time
       endif
     endif
     call obtain_interpolator_time_slices   &
                                     (terpene_interp, terpene_time)
   endif

end subroutine atmos_SOA_time_vary


!#####################################################################

subroutine atmos_SOA_endts             


      call unset_interpolator_time_flag (gas_conc_interp)
      call unset_interpolator_time_flag (isoprene_interp)
      call unset_interpolator_time_flag (terpene_interp)


end subroutine atmos_SOA_endts



!#####################################################################

!</SUBROUTINE>

!#######################################################################
!<SUBROUTINE NAME="atmos_SOA_end">
!<OVERVIEW>
!  The destructor routine for the soa module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_SOA_end
!</TEMPLATE>
 subroutine atmos_SOA_end

      call interpolator_end (gas_conc_interp)
      call interpolator_end (isoprene_interp)
      call interpolator_end (terpene_interp)
      module_is_initialized = .FALSE.

 end subroutine atmos_SOA_end

!</SUBROUTINE>
!-----------------------------------------------------------------------
      SUBROUTINE atmos_SOA_chem(pwt,temp,pfull, phalf, dt, &
                          jday,hour,minute,second,lat,lon, &
                          SOA, SOA_dt, Time,Time_next,is,ie,js,je,kbot)

! ****************************************************************************
      real, intent(in),    dimension(:,:,:)          :: pwt
      real, intent(in),    dimension(:,:,:)          :: temp,pfull,phalf
      real, intent(in)                               :: dt
      integer, intent(in)                            :: jday, hour,minute,second
      real, intent(in),  dimension(:,:)              :: lat, lon  ! [radian]
      real, intent(in),    dimension(:,:,:)          :: SOA
      real, intent(out),   dimension(:,:,:)          :: SOA_dt
      type(time_type), intent(in)                    :: Time, Time_next
      integer, intent(in),  dimension(:,:), optional :: kbot
      integer, intent(in)                            :: is,ie,js,je
! Working vectors
      real, dimension(size(SOA,1),size(SOA,2),size(SOA,3)) :: &
               SOA_chem, OH_conc, C4H10_conc
      real, dimension(size(SOA,1),size(SOA,2)) :: &
               SOA_prod, &
               xu, dayl, h, hl, hc, hred, fac_OH, fact_OH, &
	       isoprene_emis, terpene_emis

      real, parameter                            :: c4h10_SOA_yield = 0.1
      real, parameter                            :: A0 = 0.006918
      real, parameter                            :: A1 = 0.399912
      real, parameter                            :: A2 = 0.006758
      real, parameter                            :: A3 = 0.002697
      real, parameter                            :: B1 = 0.070257
      real, parameter                            :: B2 = 0.000907
      real, parameter                            :: B3 = 0.000148
      real                                       :: decl, hd, x
      integer :: i,j,k,id,jd,kd
      integer                                    :: istep, nstep
! Local grid sizes
      id=size(SOA,1); jd=size(SOA,2); kd=size(SOA,3)

      OH_conc(:,:,:)=0.  ! molec/cm3
      call interpolator(gas_conc_interp, Time, phalf, OH_conc, &
                       trim(gas_conc_name(1)), is, js)

      C4H10_conc(:,:,:)=0.0
      call interpolator(gas_conc_interp, Time, phalf, C4H10_conc, &
                       trim(gas_conc_name(2)), is, js)
      C4H10_conc(:,:,:)=C4H10_conc(:,:,:)*WTM_C4H10/WTMAIR

      x = 2. *pi *float(jday-1)/365.
      decl = A0 - A1*cos(  X) + B1*sin(  X) - A2*cos(2.*X) + B2*sin(2.*X) &
           - A3*cos(3.*X) + B3*sin(3.*X)
      xu(:,:) = -tan(lat(:,:))*tan(decl)
      where ( xu > -1 .and. xu < 1 ) dayl=acos(xu)/pi
      where ( xu <= -1 ) dayl = 1.
      where ( xu >= 1 ) dayl = 0.
!   Calculate normalization factors for OH and NO3 such that
!   the diurnal average respect the monthly input values.
      hd=0.
      fact_OH(:,:)  = 0.
      nstep = int(24.*3600./dt)
      do istep=1,nstep
        hd=hd+dt/3600./24.
        hl(:,:) = pi*(1.-dayl(:,:))
        hc(:,:) = pi*(1.+dayl(:,:))
        h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
        where ( h.ge.hl .and. h.lt.hc )
! Daytime
          hred=(h-hl)/(hc-hl)
          fact_OH  = fact_OH + amax1(0.,sin(pi*hred)/2.)/nstep
        endwhere
      enddo


      hd=amax1(0.,amin1(1.,(hour+minute/60.+second/3600.)/24.))
      hl(:,:) = pi*(1.-dayl(:,:))
      hc(:,:) = pi*(1.+dayl(:,:))
      h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
      fac_OH(:,:)  = 0.
      where ( h.ge.hl .and. h.lt.hc )
! Daytime
          hred=(h-hl)/(hc-hl)
          fac_OH  = amax1(0.,sin(pi*hred)/2.)/fact_OH
      elsewhere
! Nightime
          fac_OH  = 0.
      endwhere

!----------------------------------------------------------------------
!    SOA_dt initially contains chemical production (pseudo-emission added later)
!----------------------------------------------------------------------
      do i=1,id
        do j=1,jd
          do k=1,kd 
            SOA_dt(i,j,k) = 1.55E-11 * exp( -540./temp(i,j,k) ) * c4h10_SOA_yield &
                * C4H10_conc(i,j,k)*OH_conc(i,j,k)*fac_oh(i,j)
          enddo
        enddo
      enddo

      SOA_chem(:,:,:)=SOA_dt(:,:,:)*pwt(:,:,:)

!----------------------------------------------------------------------
!    Pseudo-emission of SOA from biogenic VOCs,
!    ... add to SOA_dt (in lowest model layer)
!----------------------------------------------------------------------
      isoprene_emis(:,:)    = 0.0
      if ( trim(isoprene_source).ne. ' ') then
        call interpolator(isoprene_interp, isoprene_time, isoprene_emis, &
                          trim(isoprene_input_name(1)), is, js)
      endif
      isoprene_emis(:,:) = isoprene_emis(:,:) * isoprene_factor * isoprene_SOA_yield

      terpene_emis(:,:)    = 0.0
      if ( trim(terpene_source).ne. ' ') then
        call interpolator(terpene_interp, terpene_time, terpene_emis, &
                          trim(terpene_input_name(1)), is, js)
      endif
      terpene_emis(:,:) = terpene_emis(:,:) * terpene_factor * terpene_SOA_yield

      SOA_dt(:,:,kd) = SOA_dt(:,:,kd) &
                     + (isoprene_emis(:,:) + terpene_emis(:,:)) / pwt(:,:,kd)

      if (id_SOA_chem > 0) then
        used = send_data ( id_SOA_chem, &
              SOA_chem, Time_next,is_in=is,js_in=js,ks_in=1)
      endif

! column production of SOA 


      SOA_prod = 0.
      do k=1,kd
        SOA_prod = SOA_prod +  SOA_chem(:,:,k)
      end do

      if (id_SOA_chem_col > 0) then
        used = send_data ( id_SOA_chem_col, &
                           SOA_prod, Time_next,is_in=is,js_in=js)
      endif

      if (id_chepsoa > 0) then
        used = send_data ( id_chepsoa, SOA_prod, Time_next, is_in=is,js_in=js)
      endif

      if (id_SOA_isoprene > 0) then
        used = send_data ( id_SOA_isoprene, &
              isoprene_emis, Time_next,is_in=is,js_in=js)
      endif
      if (id_SOA_terpene > 0) then
        used = send_data ( id_SOA_terpene, &
              terpene_emis, Time_next,is_in=is,js_in=js)
      endif
      if (id_SOA_biogenic > 0) then
        used = send_data ( id_SOA_biogenic, &
              isoprene_emis+terpene_emis, Time_next,is_in=is,js_in=js)
      endif

end subroutine atmos_SOA_chem


end module atmos_SOA_mod
