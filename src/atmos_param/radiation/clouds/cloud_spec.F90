                 module cloud_spec_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    cloud_spec_mod defines the variables that are used in a partic-
!    ular cloud parameterization to specify the cloud location, cloud
!    type and cloud magnitude for the active cloud parameterization(s).
! </OVERVIEW>
! <DESCRIPTION>
!    if microphysically-based radiative properties are desired, then
!    cloud_spec_mod also provides the microphysical parameters used in
!    determining the radiative properties, either from the cloud scheme
!    itself if they are present, or from a prescribed formula based on
!    prescribed water paths for high, middle and low clouds.
! </DESCRIPTION>

!   shared modules:

use time_manager_mod,         only: time_type, time_manager_init, &
                                    set_time, operator (+)
use mpp_mod,                  only: input_nml_file
use fms_mod,                  only: open_namelist_file, mpp_pe, &
                                    mpp_root_pe, stdlog,  fms_init, &
                                    write_version_number, file_exist, & 
                                    check_nml_error, error_mesg,   &
                                    FATAL, NOTE, close_file, stdout
use tracer_manager_mod,       only:         &
!                                   tracer_manager_init,  &
                                    get_tracer_index, NO_TRACER
use field_manager_mod,        only:       &
                                    field_manager_init, &
                                    MODEL_ATMOS
use data_override_mod,        only: data_override
use random_number_streams_mod, only: random_number_streams_init, &
                                     get_random_number_streams, &
                                     random_number_streams_end
use random_numbers_mod,    only:  randomNumberStream,   &
                                  getRandomNumbers
use constants_mod,         only : radian, RDGAS

use aerosol_types_mod,        only: aerosol_type

! atmos param modules:

use physics_radiation_exch_mod, only: clouds_from_moist_block_type, &
                                      exchange_control_type

! cloud radiation modules

use cloudrad_types_mod,       only: cld_specification_type, &
                                    microphysics_type,  &
                                    cloudrad_control_type

use strat_clouds_W_mod,       only: strat_clouds_W_init,   &
                                    strat_clouds_amt, strat_clouds_W_end

use donner_deep_clouds_W_mod, only: donner_deep_clouds_W_init, &
                                    donner_deep_clouds_amt, &
                                    donner_deep_clouds_W_end

use uw_clouds_W_mod,          only: uw_clouds_W_init, &
                                    uw_clouds_amt, &
                                    uw_clouds_W_end

!BW use rh_based_clouds_mod,      only: rh_based_clouds_init,  &
!BW                                     rh_clouds_amt, &
!BW                                     rh_based_clouds_end
                                 
!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    cloud_spec_mod defines the variables that are used in a partic-
!    ular cloud parameterization to specify the cloud location, cloud
!    type and cloud magnitude for the active cloud parameterization(s).
!    if microphysically-based radiative properties are desired, then
!    cloud_spec_mod also provides the microphysical parameters used in
!    determining the radiative properties, either from the cloud scheme
!    itself if they are present, or from a prescribed formula based on
!    prescribed water paths for high, middle and low clouds.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         cloud_spec_init, cloud_spec, cloud_spec_end

private    &

!  called from cloud_spec:
         microphys_presc_conc,  &
         combine_cloud_properties


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)  ::      &
              cloud_type_form = '     ' ! cloud parameterization being 
                                        ! used; either 'strat', 'rh', 
                                        ! 'deep',  'stratdeep', 'stratuw',
                                        ! 'stratdeepuw', 'uw', 'deepuw
                                        !  or 'none'       
real :: wtr_cld_reff=10.                ! assumed cloud drop efective
                                        ! radius [ microns ]  
real :: ice_cld_reff=50.                ! assumed ice cloud effective
                                        ! size [ microns ]
real :: rain_reff=250.                  ! assumed rain drop effective
                                        ! radius [ microns ]
character(len=16) :: overlap_type = 'random'    
                                        ! cloud overlap assumption; 
                                        ! allowable values are 'random'
                                        ! or 'max-random'  
logical :: doing_data_override=.false.
logical :: do_fu2007 = .false.
logical :: do_rain   = .false. !sjl
logical :: do_snow   = .false. !miz
logical :: do_graupel  = .false. !sjl

logical   :: do_stochastic_clouds = .false.

logical :: ignore_donner_cells = .false.! when set to .true., the effects 
                                        ! of donner cell clouds in the
                                        ! radiation code are ignored

logical :: use_cloud_tracers_in_radiation = .true.
                               ! if true, use lsc cloud tracer fields
                               ! in radiation (these transported on
                               ! current step, will have non-realizable
                               ! total cloud areas at some points); if
                               ! false, then use balanced (realizable)
                               ! fields saved at end of last step
                               ! only an issue when both lsc and conv
                               ! clouds are active (AM3)

logical :: reproduce_ulm = .true. 


namelist /cloud_spec_nml / cloud_type_form, wtr_cld_reff,   &
                           ice_cld_reff, rain_reff, overlap_type, &
                           doing_data_override, do_fu2007,    &
                           do_rain, do_snow, do_graupel, &
                           do_stochastic_clouds, &
                           ignore_donner_cells, &
                           use_cloud_tracers_in_radiation, &
                           reproduce_ulm

!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!--------------------------------------------------------------------
!    assumed water paths.
!--------------------------------------------------------------------
real   ::  lwpath_hi  = 6.313929   ! assumed water path for high clouds
                                   ! [ grams / m**2 ]
real   ::  lwpath_mid = 18.94179   ! assumed water path for middle 
                                   ! clouds [ grams / m**2 ]
real   ::  lwpath_low = 75.76714   ! assumed water path for low clouds
                                   ! [ grams / m**2 ]

!---------------------------------------------------------------------
!    logical  flags.

logical :: module_is_initialized = .false.   ! module initialized ?

!---------------------------------------------------------------------
!    time-step related constants.

integer :: num_pts       !  number of grid columns processed so far that
                         !  have cloud data present (used to identify
                         !  module coldstart condition)
integer :: tot_pts       !  total number of grid columns in the 
                         !  processor's domain

!---------------------------------------------------------------------
!     indices for cloud tracers

integer :: nql           ! tracer index for liquid water
integer :: nqi           ! tracer index for ice water
integer :: nqa           ! tracer index for cloud area
integer :: nqn           ! tracer index for cloud droplet number
integer :: nqni          ! tracer index for ice crystal number
integer :: nqr, nqs, nqg ! tracer index for rainwat, snowwat and graupel           


!----------------------------------------------------------------------
!     miscellaneous variables:

!BW integer :: num_slingo_bands  ! number of radiative bands over which 
!BW                              ! cloud optical depth is calculated in the
!BW                              ! gordon diag_cloud parameterization

integer :: id, jd, kmax

type(time_type) :: Radiation_time_step  ! saved for data override

logical :: doing_prog_clouds

!----------------------------------------------------------------------
!----------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="cloud_spec_init">
!  <OVERVIEW>
!   Contructor of cloud_spec_package module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Contructor of cloud_spec_package module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec_init ( pref, lonb, latb, axes, Time)
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!   reference pressure levels containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   the longitude array of the model grid box corners
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   the latitude array of the model grid box corners
!  </IN>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
! 
subroutine cloud_spec_init (Exch_ctrl, pref, lonb, latb, axes, Time,   &
                            rad_time_step, Cldrad_control)

!---------------------------------------------------------------------
!    cloud_spec_init is the constructor for cloud_spec_mod.
!---------------------------------------------------------------------

type(exchange_control_type), intent(inout) :: Exch_ctrl
real, dimension(:,:),        intent(in)    ::  pref        
real, dimension(:,:),        intent(in)    ::  lonb, latb
integer, dimension(4),       intent(in)    ::  axes
type(time_type),             intent(in)    ::  Time
integer,                     intent(in)    ::  rad_time_step
type(cloudrad_control_type), intent(inout) ::  Cldrad_control

!-------------------------------------------------------------------
!    intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      array of model longitudes at cell corners [ radians ]
!       latb      array of model latitudes at cell corners [radians]
!       axes      diagnostic variable axes
!       Time      current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:
 
      integer   ::   unit, ierr, io, logunit
      integer   ::   ndum, i, j, ii, jj
      

!--------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      ierr     error code
!      io       error status returned from io operation  
!      ndum     dummy argument needed for call to field_manager_init
!
!--------------------------------------------------------------------
 
!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call field_manager_init (ndum)
!  not yet compliant:
!     call tracer_manager_init  ! not public
 
!---------------------------------------------------------------------
!    read namelist.
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=cloud_spec_nml, iostat=io)
      ierr = check_nml_error(io,"cloud_spec_nml")
#else
!---------------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read (unit, nml=cloud_spec_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'cloud_spec_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
           write (logunit, nml=cloud_spec_nml)

      id = size(lonb,1) - 1
      jd = size(latb,2) - 1
      kmax = size(pref,1) - 1

!-----------------------------------------------------------------------
!    define output field.
!-----------------------------------------------------------------------
      Exch_ctrl%cloud_type_form = cloud_type_form

!--------------------------------------------------------------------
!    verify a valid type of cloud overlap. set logical variables
!    based on the namelist value.
!--------------------------------------------------------------------
      if (trim(overlap_type) == 'random') then
        Cldrad_control%do_random_overlap = .true.
      else if (trim(overlap_type) == 'max-random') then
        Cldrad_control%do_max_random_overlap = .true.
      else
        call error_mesg ('cloud_spec_mod',  &
         ' invalid specification of overlap_type', FATAL)
      endif

      doing_prog_clouds = Exch_ctrl%doing_prog_clouds

!-------------------------------------------------------------------
!    set the variables indicating that the above control variables have
!    been set.
!--------------------------------------------------------------------
!BW   Cldrad_control%do_random_overlap_iz = .true.
!BW   Cldrad_control%do_max_random_overlap_iz = .true.

!--------------------------------------------------------------------
!    save the flags indicating whether stochastic clouds are to be
!    used.
!--------------------------------------------------------------------
      Cldrad_control%do_stochastic_clouds = do_stochastic_clouds
!BW   Cldrad_control%do_stochastic_clouds_iz = .true.

!--------------------------------------------------------------------
!    if stochastic clouds is active, be sure that the
!    cloud_generator module has been initialized.
!--------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then 
          call random_number_streams_init ( lonb, latb, Cldrad_control )
      endif

!-------------------------------------------------------------------
!    verify that the nml variable cloud_type_form specifies a valid
!    cloud parameterization. set the appropriate logical control
!    variable(s) to .true.. call the constructor modules for the
!    specific cloud scheme(s) requested.
!-------------------------------------------------------------------
      if (trim(cloud_type_form) == 'strat')  then

!-------------------------------------------------------------------
!    cloud fractions, heights are predicted by the model based on klein 
!    parameterization.
!-------------------------------------------------------------------
         Cldrad_control%do_strat_clouds = .true.

!-------------------------------------------------------------------
!    cloud fractions, heights are diagnosed based on model relative 
!    humidity.
!-------------------------------------------------------------------
!BW   else if (trim(cloud_type_form)  == 'rh')   then
!BW      Cldrad_control%do_rh_clouds = .true.
!BW      call rh_based_clouds_init 

!-------------------------------------------------------------------
!    cloud fractions, heights are predicted by the donner deep cloud 
!    (cell cloud, anvil cloud) scheme.
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'deep')  then
         Cldrad_control%do_donner_deep_clouds = .true.

!------------------------------------------------------------------
!    cloud fractions, heights are provided by the uw_conv shallow
!    convection scheme  
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'uw')  then
         Cldrad_control%do_uw_clouds = .true.

!-------------------------------------------------------------------
!    cloud fractions, heights are a combination of the donner
!    deep cloud (cell cloud, anvil cloud) and klein large-scale cloud
!    parameterizations.
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'stratdeep')  then
         Cldrad_control%do_strat_clouds = .true.
         Cldrad_control%do_donner_deep_clouds = .true.

!-------------------------------------------------------------------
!    cloud fractions, heights are provided by the donner deep convection
!    (cell cloud, anvil cloud) and uw_conv shallow convection
!    cloud parameterizations.
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'deepuw')  then
         Cldrad_control%do_donner_deep_clouds = .true.
         Cldrad_control%do_uw_clouds = .true.

!-------------------------------------------------------------------
!    cloud fractions, heights are provided by the klein large-scale
!    and uw_conv shallow convection cloud parameterizations.
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'stratuw')  then
         Cldrad_control%do_strat_clouds = .true.
         Cldrad_control%do_uw_clouds = .true.

!-------------------------------------------------------------------
!    cloud fractions, heights are provided by the klein large-scale
!    the donner deep convection (cell cloud, anvil cloud) and the
!    uw_conv shallow convection cloud parameterizations.
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'stratdeepuw')  then
         Cldrad_control%do_strat_clouds = .true.
         Cldrad_control%do_donner_deep_clouds = .true.
         Cldrad_control%do_uw_clouds = .true.

!---------------------------------------------------------------
!    model is run without clouds.
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'none')  then
         Cldrad_control%do_no_clouds = .true.

!-------------------------------------------------------------------
!    failure message if none of the above options was chosen.
!-------------------------------------------------------------------
      else
         call error_mesg ('cloud_spec_mod',  &
              'invalid cloud_type_form specified', FATAL)
      endif  ! (strat)

!--------------------------------------------------------------------
!    initialize cloud schemes that are used
!--------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then
         call strat_clouds_W_init(latb, lonb, Cldrad_control, Exch_ctrl)
      endif
      if (Cldrad_control%do_donner_deep_clouds) then
         call donner_deep_clouds_W_init (pref, lonb, latb, axes, Time)
      endif
      if (Cldrad_control%do_uw_clouds) then
         call uw_clouds_W_init (Exch_ctrl)
      endif

!--------------------------------------------------------------------
!    define the dimensions of the model subdomain assigned to the 
!    processor.
!--------------------------------------------------------------------
      tot_pts = (size(latb,2)-1)*(size(lonb,1)-1)

!--------------------------------------------------------------------
!    determine if the current run is cold-starting this module. if a 
!    restart file is present, then this is not a coldstart. in that case
!    set num_pts to tot_pts so that if cloud data is not available an 
!    error message can be generated. if this is a coldstart, cloud data
!    will not be available until num_pts equals or exceeds tot_pts, so
!    continue processing without issuing an error message. 
!--------------------------------------------------------------------
      if (file_exist ('INPUT/tracer_cld_amt.res') .or.  &
          file_exist ('INPUT/strat_cloud.res') ) then
        num_pts = tot_pts
      else
        num_pts = 0
      endif

!---------------------------------------------------------------------
!    obtain the tracer indices for the strat_cloud variables when
!    running gcm.
!---------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then
          nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
          nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
          nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )

          if (do_rain) then !sjl
             nqr = get_tracer_index ( MODEL_ATMOS, 'rainwat' )
             if (nqr < 0 ) call error_mesg ('cloud_spec_mod', &
                'rainwat tracer not found, but do_rain is true', FATAL)
          end if
          if (do_snow) then !miz
             nqs = get_tracer_index ( MODEL_ATMOS, 'snowwat' )
             if (nqs < 0 ) call error_mesg ('cloud_spec_mod', &
                'snowwat tracer not found, but do_snow is true', FATAL)
          end if
          if (do_graupel) then !sjl
             nqg = get_tracer_index ( MODEL_ATMOS, 'graupel' )
             if (nqg < 0 ) call error_mesg ('cloud_spec_mod', &
                'graupel tracer not found, but do_graupel is true', FATAL)
          end if

          if (mpp_pe() == mpp_root_pe()) &
            write (logunit,'(a,3i4)') 'Stratiform cloud tracer ind&
                &ices: nql,nqi,nqa =',nql,nqi,nqa
          if (min(nql,nqi,nqa) <= 0)   &
             call error_mesg ('cloud_spec_mod', &
             'stratiform cloud tracer(s) not found', FATAL)
          if (nql == nqi .or. nqa == nqi .or. nql == nqa)   &
              call error_mesg ('cloud_spec_mod',  &
            'tracers indices cannot be the same (i.e., nql=nqi=nqa).', &
                                                              FATAL)
          nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
          if (nqn /= NO_TRACER)  then
            Cldrad_control%do_liq_num = .true.
          else
            Cldrad_control%do_liq_num = .false.
          endif
          nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )
          if (nqni /= NO_TRACER)  then
            Cldrad_control%do_ice_num = .true.
          else
            Cldrad_control%do_ice_num = .false.
          endif
      else
          Cldrad_control%do_liq_num = .false.
          Cldrad_control%do_ice_num = .false.
      endif

!BW   Cldrad_control%do_liq_num_iz = .true.
!BW   Cldrad_control%do_ice_num_iz = .true.

!---------------------------------------------------------------------
!    define the variables indicating that the cloud parameterization
!    control variables have been defined.
!---------------------------------------------------------------------
!BW   Cldrad_control%do_rh_clouds_iz = .true.
!BW   Cldrad_control%do_strat_clouds_iz = .true.
!BW   Cldrad_control%do_no_clouds_iz = .true.
!BW   Cldrad_control%do_donner_deep_clouds_iz = .true.
!BW   Cldrad_control%do_uw_clouds_iz = .true.
 
!--------------------------------------------------------------------
!    include do_fu2007 in the cloudrad_control_type variable for use
!    in other modules.
!--------------------------------------------------------------------
      Cldrad_control%using_fu2007 = do_fu2007
!BW   Cldrad_control%using_fu2007_iz = .true.     

!--------------------------------------------------------------------
!    save the radiative time step (needed when data override is active)
!--------------------------------------------------------------------
     !if (doing_data_override) then
          Radiation_time_step = set_time (rad_time_step, 0)
     !endif

!---------------------------------------------------------------------
!    mark the module initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!-------------------------------------------------------------------

end subroutine cloud_spec_init

!######################################################################
! <SUBROUTINE NAME="cloud_spec">
!  <OVERVIEW>
!    cloud_radiative_properties defines the cloud radiative properties 
!    appropriate for the radiation options that are active.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_radiative_properties defines the cloud radiative properties 
!    appropriate for the radiation options that are active.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec (is, ie, js, je, lat, z_half, z_full, Rad_time,
!                       Atmos_input, &
!                       Surface, Cld_spec, Cloud_microphys,  &
!                       Moist_clouds, r)
!  </TEMPLATE>
!  <IN NAME="is,ie,js,je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Rad_time" TYPE="time_type">
!   time at which radiation calculation is to apply
!  </IN>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </INOUT>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </INOUT>
!  <INOUT NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for 
!                        clouds associated with uw shallow convection
!  </INOUT>
!  <INOUT NAME="Surface" TYPE="Surface">
!   Surface boundary condition to radiation package
!  </INOUT>
!  <IN NAME="lsc_liquid_in" TYPE="real">
!   OPTIONAL: lsc cloud water mixing ratio
!  </IN>
!  <IN NAME="lsc_ice_in" TYPE="real">
!   OPTIONAL: cloud ice mixing ratio
!  </IN>
!  <IN NAME="lsc_area_in" TYPE="real">
!   OPTIONAL: fractional cloud area
!  </IN>
!  <IN NAME="r" TYPE="real">
!   OPTIONAL: model tracer fields on the current time step
!  </IN>
! </SUBROUTINE>
!
subroutine cloud_spec (is, ie, js, je, lat, land, z_half, z_full, Rad_time, &
                       press, pflux, temp, cloudtemp, cloudvapor, clouddeltaz, &
                       r, Cldrad_control, Cld_spec, Cloud_microphys, Aerosol, Moist_clouds_block)

!----------------------------------------------------------------------
!    cloud_spec specifies the cloud field seen by the radiation package.
!----------------------------------------------------------------------

integer,                      intent(in)             :: is, ie, js, je
real, dimension(:,:),         intent(in)             :: lat
real, dimension(:,:),         intent(in)             :: land
real, dimension(:,:,:),       intent(in)             :: z_half, z_full
type(time_type),              intent(in)             :: Rad_time
real, dimension(:,:,:),       intent(in)             :: press, pflux, &
                                                        temp, cloudtemp, &
                                                        cloudvapor, clouddeltaz
real, dimension(:,:,:,:),     intent(in)             :: r
type(cloudrad_control_type),  intent(in)             :: Cldrad_control
type(cld_specification_type), intent(inout)          :: Cld_spec    
type(microphysics_type),      intent(inout), dimension(:), allocatable :: Cloud_microphys
type(aerosol_type),           intent(in)             :: Aerosol
type(clouds_from_moist_block_type), intent(in)       :: Moist_clouds_block

!-------------------------------------------------------------------
 
!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      lat               latitude of model points  [ radians ]
!      z_half            height asl at half levels [ m ]
!      z_full            height asl at full levels [ m ]
!      Rad_time          time at which radiation calculation is to apply
!                        [ time_type (days, seconds) ] 
!
!   intent(inout) variables:
!
!      Atmos_input       atmospheric input fields on model grid,
!                        [ atmos_input_type ] 
!      Surface           variables defining the surface albedo and land
!                        fraction
!                        [ surface_type ]
!      Cld_spec          variables on the model grid which define all or
!                        some of the following, dependent on the 
!                        specific cloud parameterization: cloud optical 
!                        paths, particle sizes, cloud fractions, cloud 
!                        thickness, number of clouds in a column, 
!                        and /or cloud type (high/mid/low, ice/liq or 
!                        random/max overlap)
!                        [ cld_specification_type ]
!      Lsc_microphys     variables describing the microphysical proper-
!                        ties of the large-scale clouds
!                        [ microphysics_type ]
!      Meso_microphys    variables describing the microphysical proper-
!                        ties of the meso-scale clouds
!                        [ microphysics_type ]
!      Cell_microphys    variables describing the microphysical proper-
!                        ties of the convective cell-scale clouds
!                        [ microphysics_type ]
!
!   intent(in), optional variables:
!
!      lsc_liquid_in     cloud water mixing ratio (or specific humidity 
!                        ????) [ non-dimensional ]
!      lsc_ice_in        cloud ice mixing ratio (or specific humidity 
!                         ????) [ non-dimensional ]
!      lsc_area_in       fractional cloud area [ non-dimensional ]
!      r                 model tracer fields on the current time step
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer   :: ix, jx, kx
      logical   :: override
      logical   :: strat_data_found
      type(time_type) :: Data_time
      real, dimension (size (press,1), &
                       size (press,2), &
                       size (press,3)) :: rho
      integer   :: ncld, n
      integer   :: ncld_used
      character(len=64), dimension(size(Moist_clouds_block%Cloud_data,1)) :: scheme_names_used
      type(microphysics_type) :: Lsc_microphys

! locally define indices to make code more readable
      integer :: index_strat, index_cell, index_meso, index_shallow

!     indices for cloud schemes
!     Microphysics index (clouds actually used)

integer :: istrat, icell, imeso, ishallow

!---------------------------------------------------------------------
!   local variables:
!
!        ix      number of grid points in x direction (on processor)
!        jx      number of grid points in y direction (on processor)
!        kx      number of model layers
!        rho     atmospheric density [ kg / m**3 ]
!        ncld    number of cloud schemes
!
!---------------------------------------------------------------------

      ncld = size(Moist_clouds_block%Cloud_data,1)

!---------------------------------------------------------------------
!    check for the presence of known cloud schemes
!---------------------------------------------------------------------
      ncld_used = 0

     !-------------------
     ! stratiform clouds
     !-------------------
      istrat = 0
      strat_data_found = .false.
      if (Cldrad_control%do_strat_clouds) then
       ! maybe this condition can be allowed if tracers are to be used for large-scale
        if (Moist_clouds_block%index_strat == 0) call error_mesg ('cloud_spec_mod', &
             'stratiform cloud properties not found when &
             &stratiform clouds requested', FATAL)
        ncld_used = 1
        scheme_names_used(ncld_used) = 'strat_cloud'
        istrat = ncld_used
        strat_data_found = .true.
      endif

     !-----------------------------------------
     ! check for donner deep cloud input data
     !-----------------------------------------
      imeso = 0
      icell = 0
      if (Cldrad_control%do_donner_deep_clouds) then
        if (Moist_clouds_block%index_donner_meso == 0 .or. Moist_clouds_block%index_donner_cell == 0) & 
                   call error_mesg ('cloud_spec_mod',  &
                 'donner meso and cell properties not found when &
                 &donner clouds requested', FATAL)
        scheme_names_used(ncld_used+1) = 'donner_meso'
        scheme_names_used(ncld_used+2) = 'donner_cell'
        imeso = ncld_used+1
        icell = ncld_used+2
        ncld_used = ncld_used+2
      endif

     !-------------------------------------
     ! check for shallow cloud input data
     !-------------------------------------
      ishallow = 0
      if (Cldrad_control%do_uw_clouds) then
        if (Moist_clouds_block%index_uw_conv == 0) call error_mesg ('cloud_spec_mod',  &
                 'shallow cloud properties not found when &
                        &uw shallow clouds requested', FATAL)
        ncld_used = ncld_used+1
        scheme_names_used(ncld_used) = 'uw_conv'
        ishallow = ncld_used
      endif

      !! DEBUGGING !!
!----------------------------------------------------------------------
!    define model dimensions.
!----------------------------------------------------------------------
      ix = size(press,1)
      jx = size(press,2)
      kx = size(press,3)

!----------------------------------------------------------------------
!    allocate and initialize the arrays contained in the structures
!    used to specify the cloud amounts, types and locations and
!    the microphysical parameters.
!----------------------------------------------------------------------
      allocate(Cloud_microphys(ncld_used))
      do n = 1, ncld_used
        call Cloud_microphys(n)%alloc ( ix, jx, kx, &
                                scheme_names_used(n), Cldrad_control)
      enddo

      call Cld_spec%alloc (ix, jx, kx, Cldrad_control)

!---------------------------------------------------------------------
!    define the cloud_water, cloud_ice and cloud_area components of 
!    Cld_spec.
!---------------------------------------------------------------------
      if (Moist_clouds_block%index_strat > 0 .and. .not. use_cloud_tracers_in_radiation) then
        index_strat = Moist_clouds_block%index_strat
        Cld_spec%cloud_ice   = Moist_clouds_block%Cloud_data(index_strat)%ice_amt
        Cld_spec%cloud_water = Moist_clouds_block%Cloud_data(index_strat)%liquid_amt
        Cld_spec%cloud_area  = Moist_clouds_block%Cloud_data(index_strat)%cloud_area
        if (Cldrad_control%do_liq_num) then
            Cld_spec%cloud_droplet = Moist_clouds_block%Cloud_data(index_strat)%droplet_number
        endif
        if (Cldrad_control%do_ice_num) then
            Cld_spec%cloud_ice_num = Moist_clouds_block%Cloud_data(index_strat)%ice_number
        endif
        Cld_spec%snow       = Moist_clouds_block%Cloud_data(index_strat)%snow
        Cld_spec%rain       = Moist_clouds_block%Cloud_data(index_strat)%rain
        Cld_spec%snow_size  = Moist_clouds_block%Cloud_data(index_strat)%snow_size
        Cld_spec%rain_size  = Moist_clouds_block%Cloud_data(index_strat)%rain_size
      endif

!----------------------------------------------------------------------
!    if a cloud scheme is activated (in contrast to running without any
!    clouds), call the appropriate subroutine to define the cloud
!    location, type, amount or whatever other arrays the particular 
!    parameterization uses to specify its clouds. if the model is being
!    run with do_no_clouds = .true., exit from this routine, leaving
!    the cloud specification variables as they were initialized (to a
!    condition of no clouds).
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then

!---------------------------------------------------------------------
!    if the rh diagnostic cloud scheme is active, call rh_clouds_amt
!    to define the needed cloud specification variables.
!---------------------------------------------------------------------
!BW     if (Cldrad_control%do_rh_clouds) then
!BW       call rh_clouds_amt (is, ie, js, je, press, lat,  &
!BW                           Cld_spec)
!BW     endif ! (do_rh_clouds)

!--------------------------------------------------------------------
!    if klein prognostic clouds are active, call strat_clouds_amt to 
!    obtain the needed cloud specification variables.
!--------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds) then

!---------------------------------------------------------------------
!    if the gcm is being executed, call strat_cloud_avg to obtain the
!    appropriate (either instantaneous or time-averaged) values of
!    cloud water, cloud ice and cloud fraction. if the sa_gcm or the
!    standalone columns mode is being executed with the strat cloud
!    option, then values for the cloud water, cloud ice and when needed
!    cloud area have been input as optional arguments to this sub-
!    routine.
!---------------------------------------------------------------------
          if (Moist_clouds_block%index_strat > 0 .and. .not. use_cloud_tracers_in_radiation) then

            if (Cld_spec%cloud_area(1,1,1) == -99.) then
                Cld_spec%cloud_water(:,:,:) = r(:,:,:,nql)
                Cld_spec%cloud_ice  (:,:,:) = r(:,:,:,nqi)
                Cld_spec%cloud_area (:,:,:) = r(:,:,:,nqa)
                if (Cldrad_control%do_liq_num) then
                  Cld_spec%cloud_droplet (:,:,:) = r(:,:,:,nqn)
                endif
            endif

            if (Cld_spec%cloud_ice_num(1,1,1) == -99.) then
                if (Cldrad_control%do_ice_num) then
                  Cld_spec%cloud_ice_num (:,:,:) = r(:,:,:,nqni)
                endif
            endif

          else  ! (present (lsc_liquid_in))

              Cld_spec%cloud_water(:,:,:) = r(:,:,:,nql)
              Cld_spec%cloud_ice  (:,:,:) = r(:,:,:,nqi)
              Cld_spec%cloud_area (:,:,:) = r(:,:,:,nqa)
              if (Cldrad_control%do_liq_num) then
                Cld_spec%cloud_droplet (:,:,:) = r(:,:,:,nqn)
              endif
              if (Cldrad_control%do_ice_num) then
                Cld_spec%cloud_ice_num (:,:,:) = r(:,:,:,nqni)
              endif
          endif ! (present(lsc_liquid_in))


            if (do_rain) then !sjl
              Cld_spec%cloud_water(:,:,:) = Cld_spec%cloud_water(:,:,:)+r(:,:,:,nqr)
            end if
            if (do_snow) then !miz
              Cld_spec%cloud_ice(:,:,:) = Cld_spec%cloud_ice(:,:,:)+r(:,:,:,nqs)
            end if
            if (do_graupel) then !SJL
              Cld_spec%cloud_ice(:,:,:) = Cld_spec%cloud_ice(:,:,:)+r(:,:,:,nqg)
            end if

!---------------------------------------------------------------------
!    if the cloud input data is to be overriden, define the time slice
!    of data which is to be used. allocate storage for the cloud data.
!---------------------------------------------------------------------
          if (doing_data_override) then
            Data_time = Rad_time + Radiation_time_step
 
!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    water data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current 
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qlnew', Cld_spec%cloud_water,   &
                                Data_time, override=override,           &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'ql => cloud_water not overridden successfully', FATAL)
            endif

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    ice data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current 
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qinew', Cld_spec%cloud_ice,   &
                                Data_time, override=override,         &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'qi => cloud_ice   not overridden successfully', FATAL)
            endif

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    fraction data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qanew', Cld_spec%cloud_area,   &
                                Data_time, override=override,         &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
               'qa => cloud_area not overridden successfully', FATAL)
            endif
            strat_data_found = .true.
          endif ! (doing_override)

!---------------------------------------------------------------------
!    if values for the cloud variables have been successfully obtained,
!    call strat_clouds_amt to define the appropriate cloud specification
!    variables.
!---------------------------------------------------------------------
          if (strat_data_found) then
            call strat_clouds_amt (is, ie, js, je, Rad_time, &
                                   pflux, press, cloudtemp, &
                                   cloudvapor(:,:,:)/(1.0+cloudvapor(:,:,:)), &
                                   land, Cldrad_control, &
                                   Cld_spec, Cloud_microphys(istrat), Aerosol)

!----------------------------------------------------------------------
!    cloud data was not successfully obtained.
!    if this is not the coldstart step, write an error message and 
!    stop execution.
!----------------------------------------------------------------------
          else 
            if (num_pts >= tot_pts) then
              call error_mesg ('cloud_spec_mod',  &
                     'no strat cloud data available', FATAL)

!----------------------------------------------------------------------
!    if this is the coldstart step, retain the input values corres-
!    ponding to no clouds, increment the points counter, and continue. 
!----------------------------------------------------------------------
            else
!$OMP ATOMIC UPDATE
              num_pts = num_pts + size(press,1)*size(press,2)
            endif
          endif

        else ! (do_strat_clouds)
!----------------------------------------------------------------------
!    only allocate part of the variables in this derived type
!    when strat_clouds is not active
!       (WHERE DOES THIS GET DEALLOCATED?)
!----------------------------------------------------------------------
          allocate(Lsc_microphys%conc_drop(id,jd,kmax))
          allocate(Lsc_microphys%conc_ice (id,jd,kmax))
          allocate(Lsc_microphys%size_drop(id,jd,kmax))
          allocate(Lsc_microphys%size_ice (id,jd,kmax))
          allocate(Lsc_microphys%size_rain(id,jd,kmax))
          Lsc_microphys%conc_drop = 0.0
          Lsc_microphys%conc_ice  = 0.0
          Lsc_microphys%size_drop = 0.0
          Lsc_microphys%size_ice  = 0.0
          Lsc_microphys%size_rain = 0.0

        endif ! (do_strat_clouds)

!--------------------------------------------------------------------
!    since donner_deep_clouds may be active along with strat clouds, 
!    the associated properties are determined outside of the above loop.
!    these properties are placed in Cell_microphys and Meso_microphys.
!----------------------------------------------------------------------
        if (Cldrad_control%do_donner_deep_clouds) then
          index_cell = Moist_clouds_block%index_donner_cell
          index_meso = Moist_clouds_block%index_donner_meso

          call donner_deep_clouds_amt (is, ie, js, je,  &
                           Moist_clouds_block%Cloud_data(index_cell)%cloud_area,  &
                           Moist_clouds_block%Cloud_data(index_cell)%liquid_amt,  &
                           Moist_clouds_block%Cloud_data(index_cell)%liquid_size, &
                           Moist_clouds_block%Cloud_data(index_cell)%ice_amt,     &
                           Moist_clouds_block%Cloud_data(index_cell)%ice_size,    &
                           Moist_clouds_block%Cloud_data(index_cell)%droplet_number, &
                           Moist_clouds_block%Cloud_data(index_meso)%cloud_area,  &
                           Moist_clouds_block%Cloud_data(index_meso)%liquid_amt,  &
                           Moist_clouds_block%Cloud_data(index_meso)%liquid_size, &
                           Moist_clouds_block%Cloud_data(index_meso)%ice_amt,     &
                           Moist_clouds_block%Cloud_data(index_meso)%ice_size,    &
                           Moist_clouds_block%Cloud_data(index_meso)%droplet_number, &
                           Moist_clouds_block%Cloud_data(index_meso)%nsum_out,       &
                           Cloud_microphys(icell), Cloud_microphys(imeso)  )

!---------------------------------------------------------------------
!    convert the cloud and ice amounts from kg(h2o) / kg(air) to 
!    g(h2o) / m**3, as required for use in the microphys_rad routines
!    which compute cloud radiative properties.
!---------------------------------------------------------------------
          rho(:,:,:) = press(:,:,1:kx)/(RDGAS*temp(:,:,1:kx))
          Cloud_microphys(icell)%conc_drop = 1.0e03*rho*Cloud_microphys(icell)%conc_drop
          Cloud_microphys(icell)%conc_ice  = 1.0e03*rho*Cloud_microphys(icell)%conc_ice
          Cloud_microphys(imeso)%conc_drop = 1.0e03*rho*Cloud_microphys(imeso)%conc_drop
          Cloud_microphys(imeso)%conc_ice  = 1.0e03*rho*Cloud_microphys(imeso)%conc_ice
        endif

!--------------------------------------------------------------------
!    since uw_clouds may be active along with strat clouds and / or 
!    donner deep clouds, the associated properties are determined 
!    outside of the above loop. these properties are placed in  
!    Shallow_microphys.
!----------------------------------------------------------------------
        if (Cldrad_control%do_uw_clouds) then
          index_shallow = Moist_clouds_block%index_uw_conv
          call uw_clouds_amt (is, ie, js, je,  &
                           Moist_clouds_block%Cloud_data(index_shallow)%cloud_area,     &
                           Moist_clouds_block%Cloud_data(index_shallow)%liquid_amt,     &
                           Moist_clouds_block%Cloud_data(index_shallow)%ice_amt,        &
                           Moist_clouds_block%Cloud_data(index_shallow)%droplet_number, &
                           Moist_clouds_block%Cloud_data(index_shallow)%ice_number,     &
                           land, press, cloudtemp, Cldrad_control, Cloud_microphys(ishallow) )
        endif

!---------------------------------------------------------------------
!    obtain the microphysical properties (sizes and concentrations) if
!    a prescribed microphysics scheme is active. 
!---------------------------------------------------------------------
        if (Cldrad_control%do_presc_cld_microphys) then
          call microphys_presc_conc (is, ie, js, je,   &
                                     clouddeltaz, cloudtemp, &
                                     Cld_spec, Lsc_microphys)
        endif

!---------------------------------------------------------------------
!    call combine_cloud_properties to combine (if necessary) the cloud 
!    properties from multiple cloud types (large-scale, donner deep,
!    uw shallow) into a single set for use by the radiation package. 
!    this is only needed when microphysically-based properties are 
!    present, and when either strat clouds, donner deep and / or uw
!    shallow clouds is activated.
!---------------------------------------------------------------------
!BW     if (Cldrad_control%do_sw_micro .or. Cldrad_control%do_lw_micro) then
          if (ncld_used > 0) then
              call combine_cloud_properties ( is, js,  &
                                             temp(:,:,1), Rad_time, &
                                             Cldrad_control, Cloud_microphys, Cld_spec)
          endif
!BW     endif


!--------------------------------------------------------------------
!    if microphysics is active and strat_clouds is not, define the water
!    paths (in units of kg / m**2).  if strat_clouds is active, these 
!    values will have already been defined. when microphysics is active,
!    define the effective sizes for the liquid and ice particles.
!--------------------------------------------------------------------
!BW   if (Cldrad_control%do_lw_micro .or.    &
!BW       Cldrad_control%do_sw_micro)  then
    !BW if (.not. Cldrad_control%do_strat_clouds) then

        if (Cldrad_control%do_presc_cld_microphys) then
          Cld_spec%lwp = 1.0E-03*Lsc_microphys%conc_drop(:,:,:)*clouddeltaz(:,:,:)
          Cld_spec%iwp = 1.0E-03*Lsc_microphys%conc_ice(:,:,:)*clouddeltaz(:,:,:)
          Cld_spec%reff_liq_micro = Lsc_microphys%size_drop
          Cld_spec%reff_ice_micro = Lsc_microphys%size_ice
          deallocate(Lsc_microphys%conc_drop, &
                     Lsc_microphys%conc_ice,  &
                     Lsc_microphys%size_drop, &
                     Lsc_microphys%size_ice,  &
                     Lsc_microphys%size_rain)
        endif

        if (Cldrad_control%do_strat_clouds) then
          Cld_spec%reff_liq_micro = Cloud_microphys(istrat)%size_drop
          Cld_spec%reff_ice_micro = Cloud_microphys(istrat)%size_ice
        endif

      endif  !  (.not. do_no_clouds)
!---------------------------------------------------------------------


end subroutine cloud_spec    


!######################################################################

subroutine cloud_spec_end (Cldrad_control)

type(cloudrad_control_type), intent(in) :: Cldrad_control

!---------------------------------------------------------------------
!    cloud_spec_end is the destructor for cloud_spec_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    close the modules that were initialized by this module.
!--------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then

!-------------------------------------------------------------------
!    rh-based diagnostic clouds.
!-------------------------------------------------------------------
!BW     if (Cldrad_control%do_rh_clouds) then
!BW       call rh_based_clouds_end 

!------------------------------------------------------------------
!    cloud types which may coexist must be processed outside of if loop
!------------------------------------------------------------------
!BW     else

          if (Cldrad_control%do_strat_clouds) then
            call strat_clouds_W_end (Cldrad_control)
          endif
          if (Cldrad_control%do_donner_deep_clouds) then
            call donner_deep_clouds_W_end
          endif
          if (Cldrad_control%do_uw_clouds) then
            call uw_clouds_W_end
          endif
!BW     endif

        if (Cldrad_control%do_stochastic_clouds) then
          call random_number_streams_end
        endif

      endif  ! (not do_no_clouds)

!--------------------------------------------------------------------
!    mark the module as no longer initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------

end subroutine cloud_spec_end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!###################################################################
! <SUBROUTINE NAME="combine_cloud_properties">
!  <OVERVIEW>
!    combine_cloud_properties produces cloud specification property 
!    arrays for the total cloud field in each grid box, using as input 
!    the specification of the component cloud types that may be present
!    (large-scale, mesoscale and cell-scale).
!  </OVERVIEW>
!  <DESCRIPTION>
!    combine_cloud_properties produces cloud specification property 
!    arrays for the total cloud field in each grid box, using as input 
!    the specification of the component cloud types that may be present
!    (large-scale, mesoscale and cell-scale).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call combine_cloud_properties (Lsc_microphys, Meso_microphys,  &
!                                     Cell_microphys, Cld_spec)
!  </TEMPLATE>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </IN>
!  <IN NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </IN>
!  </IN>
!  <IN NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for 
!                        clouds associated with uw shallow convection
!  </IN>
! </SUBROUTINE>
! 
subroutine combine_cloud_properties (is, js, temp, Rad_time, &
                                     Cldrad_control, Cloud_microphys, Cld_spec)

!----------------------------------------------------------------------
!    combine_cloud_properties produces cloud specification property 
!    arrays for the total cloud field in each grid box, using as input 
!    the specification of the component cloud types that may be present
!    (large-scale, donner mesoscale and cell-scale, uw shallow).
!----------------------------------------------------------------------

integer, intent(in)  :: is, js
real, dimension(:,:), intent(in) :: temp
type(time_type), intent(in) :: Rad_time
type(cloudrad_control_type),            intent(in)    :: Cldrad_control
type(microphysics_type),  dimension(:), intent(in)    :: Cloud_microphys
type(cld_specification_type),           intent(inout) :: Cld_spec

!----------------------------------------------------------------------
!   intent(in) variables:
!
!       Lsc_microphys  microphysical specification for large-scale 
!                      clouds
!                      [ microphysics_type ]
!       Meso_microphys microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!                      [ microphysics_type ]
!       Cell_microphys microphysical specification for convective cell
!                      clouds associated with donner convection
!                      [ microphysics_type ]
!       Shallow_microphys 
!                      microphysical specification for 
!                      clouds associated with uw shallow convection
!                      [ microphysics_type ]
!
!    intent(inout) variables:
!
!       Cld_spec       variables on the model grid which define all or
!                      some of the following, dependent on the specific
!                      cloud parameterization: cloud optical paths, 
!                      particle sizes, cloud fractions, cloud thickness,
!                      number of clouds in a column, and /or cloud type 
!                      (high/mid/low, ice/liq or random/max overlap)
!                      [ cld_specification_type ]
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!    variables for folding Donner cloud properties into stochastic
!    cloud arrays
!------------------------------------------------------------------
      type(randomNumberStream),   &
                    dimension(size(Cld_spec%camtsw,1),   &
                              size(Cld_spec%camtsw,2)) :: streams
      real, &            
                    dimension(size(Cld_spec%camtsw,1),   &
                              size(Cld_spec%camtsw,2),   &       
                              size(Cld_spec%camtsw,3),   &       
                              Cldrad_control%num_lw_cloud_bands+ &
                              Cldrad_control%num_sw_cloud_bands) :: &
                                                     randomNumbers
      integer :: nn, nsubcols

      integer :: i, j, k, n, ncld
      integer :: meso, cell
      integer :: istrat, icell, imeso, ishallow

!---------------------------------------------------------------------
!    total-cloud specification properties need be defined only when
!    strat_cloud, donner_deep and/or uw shallow clouds are active.
!---------------------------------------------------------------------

      ncld = size(Cloud_microphys,1)

  ! indices for cloud types in microphysics type
      istrat=0
      icell=0
      imeso=0
      ishallow=0
      do n = 1, ncld
         if (trim(Cloud_microphys(n)%scheme_name) == 'strat_cloud') istrat = n
         if (trim(Cloud_microphys(n)%scheme_name) == 'donner_meso') imeso  = n
         if (trim(Cloud_microphys(n)%scheme_name) == 'donner_cell') icell  = n
         if (trim(Cloud_microphys(n)%scheme_name) == 'uw_conv')     ishallow = n
      enddo

!----------------------------------------------------------------------
!    define the random overlap cloud fraction as the sum of the
!    fractions of all cloud schemes
!---------------------------------------------------------------------

      if (reproduce_ulm .and. ncld == 4) then
         Cld_spec%crndlw = Cloud_microphys(istrat)%cldamt + &
                           Cloud_microphys(icell)%cldamt + &
                           Cloud_microphys(imeso)%cldamt + &
                           Cloud_microphys(ishallow)%cldamt

      else if (reproduce_ulm .and. ncld == 3) then
         Cld_spec%crndlw = Cloud_microphys(1)%cldamt + &
                           Cloud_microphys(2)%cldamt + &
                           Cloud_microphys(3)%cldamt

      else if (reproduce_ulm .and. ncld == 2) then
         Cld_spec%crndlw = Cloud_microphys(1)%cldamt + &
                           Cloud_microphys(2)%cldamt

      else if (reproduce_ulm .and. ncld == 1) then
         Cld_spec%crndlw = Cloud_microphys(1)%cldamt

      else
       ! general case: does not reproduce ulm version
         Cld_spec%crndlw = 0.0
         do n = 1, ncld
            Cld_spec%crndlw = Cld_spec%crndlw + Cloud_microphys(n)%cldamt
         enddo
      endif

!---------------------------------------------------------------------
!    randomly-overlapped clouds are being assumed for donner_deep and 
!    strat cloud module clouds. set the max overlap cloud fraction to 
!    zero, be certain that the random overlap fraction is .le. 1. after
!    the summing of the component cloud fractions, and define the total
!    cloud fraction to be used by the sw code.
!---------------------------------------------------------------------
      Cld_spec%cmxolw = 0.0
      Cld_spec%crndlw = MIN (Cld_spec%crndlw, 1.00)
      Cld_spec%camtsw = Cld_spec%crndlw

!--------------------------------------------------------------------
!    if stochastic clouds are being used, define the cloud type to be 
!    seen by the radiation code in each stochastic subcolumn.
!--------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then
        nsubcols = Cldrad_control%num_sw_cloud_bands + &
                   Cldrad_control%num_lw_cloud_bands

!--------------------------------------------------------------------
!   assign either a 1 or a 0 to each subcolumn indicating whether
!   lsc cloud is present or not. 
!--------------------------------------------------------------------
        if (istrat > 0) then
          do n=1,nsubcols
            if ( n > Cldrad_control%num_sw_cloud_bands) then
              nn = n - Cldrad_control%num_sw_cloud_bands    
            else
              nn = n + Cldrad_control%num_lw_cloud_bands    
            endif
            do k=1,size(Cld_spec%camtsw,3) ! Levels
               do j=1,size(Cld_spec%camtsw,2) ! Lons
                do i=1,size(Cld_spec%camtsw,1) ! Lats
                  if (Cloud_microphys(istrat)%stoch_cldamt(i,j,k,nn) > 0.) then
                    !----------------------------------------------
                    ! fill it in with the large-scale cloud values
                    !-----------------------------------------------
                     Cld_spec%stoch_cloud_type(i,j,k,n) = istrat
                    !Cld_spec%stoch_cloud_type(i,j,k,n) = 1
                  else
                     Cld_spec%stoch_cloud_type(i,j,k,n) = 0  
                  endif
                enddo
              enddo
            enddo
          enddo
        endif ! strat_cloud

!----------------------------------------------------------------------
!    compare the cell and meso-scale cloud amounts to a random number, 
!    and replace the large-scale cloud and clear sky assignment in each 
!    subcolumn with an assignment of cell or meso-scale clouds when the 
!    number is less than the cloud fraction. use the maximum overlap 
!    assumption. treat the random number as the location with the PDF 
!    of total water. cells are at the top of the PDF; then meso-scale 
!    anvils, then large-scale clouds and clear sky.
!------------------------------------------------------------
        if (Cldrad_control%do_donner_deep_clouds) then     
          call get_random_number_streams (is, js, Rad_time, temp, streams, perm=1)

!----------------------------------------------------------------------
!    get the random numbers to do both sw and lw at oncer.
!----------------------------------------------------------------------
          do j=1,size(Cld_spec%camtsw,2) ! Lons
            do i=1,size(Cld_spec%camtsw,1) ! Lats
              call getRandomNumbers (streams(i,j), randomNumbers(i,j,1,:))
            end do
          end do
 
!----------------------------------------------------------------------
!    here is maximum overlap. we use a 3D arrary for the random numbers
!    for flexibility.
!----------------------------------------------------------------------
          do k=2,size(Cld_spec%camtsw,3)
            randomNumbers(:,:,k,:) = randomNumbers(:,:,1,:)
          end do
 
!----------------------------------------------------------------------
!    assign cloud types, band by band
!----------------------------------------------------------------------
          if (ignore_donner_cells) then
            do n=1,nsubcols    
              where( randomNumbers(:,:,:,n) > &
                     (1. - Cloud_microphys(imeso)%cldamt)) &
                   ! assign meso-scale cloud
                     Cld_spec%stoch_cloud_type(:,:,:,n) = imeso
                   ! Cld_spec%stoch_cloud_type(:,:,:,n) = 2
            enddo

          else ! ignore_donner_cells
            do n=1,nsubcols    
               where( randomNumbers(:,:,:,n) > &
                      (1. - Cloud_microphys(icell)%cldamt - &
                            Cloud_microphys(imeso)%cldamt)) &
                   ! assign meso-scale cloud
                     Cld_spec%stoch_cloud_type(:,:,:,n) = imeso
                    !Cld_spec%stoch_cloud_type(:,:,:,n) = 2  
               where( randomNumbers(:,:,:,n) > &
                      (1. - Cloud_microphys(icell)%cldamt)) &
                   ! assign meso-scale cloud
                     Cld_spec%stoch_cloud_type(:,:,:,n) = icell
                    !Cld_spec%stoch_cloud_type(:,:,:,n) = 3  
            enddo ! n
          endif ! ignore_donner_cells
        endif ! do_donner_deep_clouds

!----------------------------------------------------------------------
!    compare the uw shallow cloud amount to a random number, and replace
!    the donner cloud, large-scale cloud or clear sky previously 
!    assigned in each subcolumn with an assignment of uw shallow cloud 
!    when the number is less than the cloud fraction. use the maximum 
!    overlap assumption. treat the random number as the location with 
!    the PDF of total water. uw shallow clouds are at the top of this 
!    PDF, then large-scale clouds and clear sky.
!------------------------------------------------------------
        if (Cldrad_control%do_uw_clouds) then     
          call get_random_number_streams (is, js, Rad_time, temp, streams, perm=2)

!----------------------------------------------------------------------
!    get the random numbers to do both sw and lw at oncer.
!----------------------------------------------------------------------
          do j=1,size(Cld_spec%camtsw,2) ! Lons
            do i=1,size(Cld_spec%camtsw,1) ! Lats
              call getRandomNumbers (streams(i,j), randomNumbers(i,j,1,:))
            end do
          end do
 
!----------------------------------------------------------------------
!    here is maximum overlap. we use a 3D arrary for the random numbers
!    for flexibility.
!----------------------------------------------------------------------
          do k=2,size(Cld_spec%camtsw,3)
            randomNumbers(:,:,k,:) = randomNumbers(:,:,1,:)
          end do
 
!----------------------------------------------------------------------
!    assign cloud type, band by band
!----------------------------------------------------------------------
          do n=1,nsubcols
            where( randomNumbers(:,:,:,n) > &
                   (1. - Cloud_microphys(ishallow)%cldamt)) &
                 ! assign uw shallow cloud
                   Cld_spec%stoch_cloud_type(:,:,:,n) = ishallow
                 ! Cld_spec%stoch_cloud_type(:,:,:,n) = 4
          enddo
        endif

!---------------------------------------------------------------------
!     define the cloud amount in each stochastic subcolumn to be either
!     1.0 if cloud is present, or 0.0 if no cloud exists.
!---------------------------------------------------------------------
        do n=1,Cldrad_control%num_sw_cloud_bands
          do k=1,size(Cld_spec%camtsw,3) ! Levels
            do j=1,size(Cld_spec%camtsw,2) ! Lons
              do i=1,size(Cld_spec%camtsw,1) ! Lats
                if (Cld_spec%stoch_cloud_type(i,j,k,n) /= 0) then  
                  Cld_spec%camtsw_band(i,j,k,n) = 1.0
                else
                  Cld_spec%camtsw_band(i,j,k,n) = 0.0
                endif
              end do
            end do
          end do
        end do
        
        do n=1,Cldrad_control%num_lw_cloud_bands
          nn = Cldrad_control%num_sw_cloud_bands + n
          do k=1,size(Cld_spec%camtsw,3) ! Levels
            do j=1,size(Cld_spec%camtsw,2) ! Lons
              do i=1,size(Cld_spec%camtsw,1) ! Lats
                if (Cld_spec%stoch_cloud_type(i,j,k,nn) /= 0) then  
                  Cld_spec%crndlw_band(i,j,k,n) = 1.0
                else
                  Cld_spec%crndlw_band(i,j,k,n) = 0.0
                endif
              end do
            end do
          end do
        end do
      endif  ! (do_stochastic)


!-------------------------------------------------------------------


end subroutine combine_cloud_properties 



!###################################################################
! <SUBROUTINE NAME="microphs_presc_conc">
!  <OVERVIEW>
!   Subroutine to determine water droplet and ice crystal based on
!   prescribed microphysics model.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine uses prescribed microphysics model to determine
!   concentrations of water droplets and ice crystals. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call microphys_presc_conc (is, ie, js, je, deltaz, temp,      &
!                                 Cld_spec, Lsc_microphys)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   starting indice of the x dimension in the physics domain
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   ending indice of the x dimension in the physics domain
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   starting indice of the y dimension in the physics domain
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   ending indice of the y dimension in the physics domain 
!  </IN>
!  <IN NAME="deltaz" TYPE="real">
!   Height of each pressure layers.
!  </IN>
!  <IN NAME="temp" TYPE="real">
!   Temperatures of pressure levels
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
! </SUBROUTINE>
!
subroutine microphys_presc_conc (is, ie, js, je, deltaz, temp,      &
                                 Cld_spec, Lsc_microphys)

!---------------------------------------------------------------------
!    microphys_presc_conc defines microphysical properties based on the
!    assumption of specified total water paths for high, middle and low 
!    clouds.
!---------------------------------------------------------------------

integer,                      intent(in)     :: is, ie, js, je
real, dimension(:,:,:),       intent(in)     :: deltaz, temp  
type(cld_specification_type), intent(in)     :: Cld_spec
type(microphysics_type),      intent(inout)  :: Lsc_microphys

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      deltaz         model vertical grid separation that is to be used
!                     for cloud calculations
!                     [meters]
!      temp           temperature at model levels (1:nlev) that is to
!                     be used in cloud calculations
!                     [ deg K ]
!      Cld_spec       cld_specification_type structure, contains var-
!                     iables defining the cloud distribution
!
!   intent(inout) variables:
!
!      Lsc_microphys  microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     large-scale clouds
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
! local variables:                                                  
!---------------------------------------------------------------------

      real,    dimension(size(temp,1), size(temp,2), size(temp,3)) :: &
                                                       conc

      integer, dimension(size(temp,1), size(temp,2)) :: &
                                                       nhi_clouds, &
                                                       nmid_clouds, &
                                                       nlow_clouds

      integer  :: i,j,k

!--------------------------------------------------------------------
!  local variables:
!
!      conc             droplet concentration  [ g / m**3 ]
!      nhi_clouds       number of layers with high clouds
!      nmid_clouds      number of layers with middle clouds
!      nlow_clouds      number of layers with low clouds
!      i,j,k            do-loop indices
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!! RSH NOTE:
!
!    THE FOLLOWING treatment of diag_cloud_mod is here as an INITIAL 
!    IMPLEMENTATION to allow compilation and model execution, and 
!    provide "reasonable ?? " values.
! 
!    Code developed but NOT YET ADDED HERE reflects a later approach. 
!    That code is available under the fez release, and will be added to
!    the repository when upgrades to the cloud-radiation modules are 
!    completed.
!
!    obtain drop and ice size and concentrations, consistent with 
!    the diag_cloud scheme. As a test case, the following is a simple 
!    specification of constant concentration and size in all boxes 
!    defined as cloudy, attempting to come close to the prescribed 
!    values used for other cloud schemes. assume ice cld thickness 
!    = 2.0 km; then conc_ice=10.0E-03 => iwp = 20 g/m^2, similar to that
!    prescribed in microphys_presc_conc. assume water cld thickness 
!    = 3.5 km; then conc_drop = 20E-03 => lwp = 70 g / m^2, similar to 
!    that prescribed in microphys_presc_conc.  use sizes as used in 
!    microphys_presc_conc (50 and 20 microns). when done, radiative 
!    boundary fluxes are "similar" to non-microphysical results
!    for test case done here, and shows reasonable sensitivity to
!    variations in concentrations.
!    AGAIN, THIS IS AN INITIAL IMPLEMENTATION FOR TESTING ONLY !!!!
!
!!! RSH
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!---------------------------------------------------------------------
!    for the non-diag_cloud_mod cases, assume that the water path is 
!    preset at fixed values (lwpath_hi, _mid, _low) for "high", "mid", 
!    "low" clouds. the lwpath in each cloud layer within "hi", "mid" 
!    "low" pressure intervals is that lwpath_... divided by the number 
!    of clouds present in that pressure interval.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    define the number of high, middle, low clouds according to
!    Wetherald's criterion.
!----------------------------------------------------------------------
        do j=1,size(Cld_spec%camtsw,2)
          do i=1,size(Cld_spec%camtsw,1)
            nhi_clouds(i,j)  = 0
            nmid_clouds(i,j) = 0
            nlow_clouds(i,j) = 0
            do k=1,size(Cld_spec%camtsw,3)
              if (Cld_spec%hi_cloud(i,j,k)) &
                               nhi_clouds(i,j)  =  nhi_clouds(i,j)  + 1
              if (Cld_spec%mid_cloud(i,j,k)) &
                               nmid_clouds(i,j) =  nmid_clouds(i,j) + 1
              if (Cld_spec%low_cloud(i,j,k))  &
                               nlow_clouds(i,j) =  nlow_clouds(i,j) + 1
            end do
          end do
        end do

!----------------------------------------------------------------------
!    compute the water substance concentration in each layer 
!    (as water path / layer geometric path).
!----------------------------------------------------------------------
        conc(:,:,:) = 0.0E+00
        do j=1,size(Cld_spec%camtsw,2)
          do i=1,size(Cld_spec%camtsw,1)
            do k=1,size(Cld_spec%camtsw,3)
              if (Cld_spec%hi_cloud(i,j,k)) then
                conc(i,j,k) = lwpath_hi/   &
                              (nhi_clouds(i,j)*deltaz(i,j,k))
              endif
              if (Cld_spec%mid_cloud(i,j,k)) then
                conc(i,j,k) = lwpath_mid/    &
                              (nmid_clouds(i,j)*deltaz(i,j,k))
              endif
              if (Cld_spec%low_cloud(i,j,k)) then
                conc(i,j,k) = lwpath_low    /                   &
                              (nlow_clouds(i,j)*deltaz(i,j,k))
              endif
            end do
          end do
        end do

!----------------------------------------------------------------------
!    split conc into conc_ice and conc_drop, depending on temperature
!    criterion (T < 273.16). assume that rain and / or snow are not
!    present.
!----------------------------------------------------------------------
        do k=1,size(Cld_spec%camtsw,3)
          do j=1,size(Cld_spec%camtsw,2)
            do i=1,size(Cld_spec%camtsw,1)
              if (temp(i,j,k) .LT. 273.16) then
                Lsc_microphys%conc_ice(i,j,k) = conc(i,j,k)
              else
                Lsc_microphys%conc_drop(i,j,k) = conc(i,j,k)
              endif
            end do
          end do
        end do

!----------------------------------------------------------------------
!    define sizes of microphysical species, using namelist values. note
!    that namelist drop and rain sizes are radii, so multiply by 2 to 
!    produce diameter, as desired for the %size_ arrays.
!----------------------------------------------------------------------
        Lsc_microphys%size_drop(:,:,:) = 2.0*wtr_cld_reff
        Lsc_microphys%size_rain(:,:,:) = 2.0*rain_reff
        Lsc_microphys%size_ice (:,:,:) = ice_cld_reff

!--------------------------------------------------------------------


end subroutine microphys_presc_conc


!#################################################################


                end module cloud_spec_mod


