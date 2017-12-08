
       module cloudrad_driver_mod

use mpp_mod,               only: input_nml_file
use fms_mod,               only: fms_init, mpp_clock_id, &
                                 mpp_clock_begin, mpp_clock_end, &
                                 CLOCK_MODULE_DRIVER, CLOCK_MODULE, &
                                 field_exist, field_size, &
                                 mpp_pe, mpp_root_pe, &
                                 open_namelist_file, stdlog, stdout, &
                                 file_exist, FATAL, WARNING, NOTE, &
                                 close_file, read_data, write_data, &
                                 write_version_number, check_nml_error,&
                                 error_mesg, mpp_chksum, &
                                 read_data, mpp_error
use time_manager_mod,      only: time_type

use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: tracer_manager_init, &
                                 get_number_tracers

use aerosol_types_mod,     only: aerosol_type

use physics_radiation_exch_mod, only: clouds_from_moist_block_type, &
                                      exchange_control_type

!  cloud radiation modules:

use cloudrad_types_mod,    only: cld_specification_type, &
                                 microphysics_type, &
                                 cloudrad_control_type

use cloud_spec_mod,        only: cloud_spec_init, &
                                 cloud_spec_end, &
                                 cloud_spec

use cloudrad_package_mod,  only: cloudrad_package_init, &
                                 cloud_radiative_properties, &
                                 cloudrad_package_end

!--------------------------------------------------------------------

implicit none 
private 

!----------------------------------------------------------------------
!------------ version number for this module --------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

!---------------------------------------------------------------------
!------ interfaces -----

public    cloudrad_driver_init, &
          cloudrad_driver, &
          cloudrad_driver_end

!--------------------------------------------------------------------
!------- namelist ---------

integer :: dummy = 0

namelist /cloudrad_driver_nml/ dummy

!---------------------------------------------------------------------
!    miscellaneous control variables:

logical ::  module_is_initialized = .false. ! module initialized?

!---------------------------------------------------------------------

                         contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine cloudrad_driver_init (Time, rad_time_step,  &
                                 lonb, latb, axes, pref, &
                                 num_sw_cloud_bands, &
                                 num_lw_cloud_bands, &
                                 Cldrad_control, Exch_ctrl )

!---------------------------------------------------------------------
!   cloudrad_driver_init is the constructor
!   for cloudrad_driver_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
type(time_type),             intent(in)    :: Time
integer,                     intent(in)    :: rad_time_step
real, dimension(:,:),        intent(in)    :: lonb, latb
integer, dimension(4),       intent(in)    :: axes
real, dimension(:,:),        intent(in)    :: pref
integer,                     intent(in)    :: num_sw_cloud_bands, &
                                              num_lw_cloud_bands
type(cloudrad_control_type), intent(inout) :: Cldrad_control
type(exchange_control_type), intent(inout) :: Exch_ctrl
!----------------------------------------------------------------------
!   intent(in) variables:
!
!       Time           current time [time_type(days, seconds)]
!       lonb           2d array of model longitudes on cell corners 
!                      [ radians ]
!       latb           2d array of model latitudes at cell corners 
!                      [ radians ]
!       axes           diagnostic variable axes
!       pref           array containing two reference pressure profiles 
!                      for use in defining transmission functions
!                      [ pascals ]
!----------------------------------------------------------------------
!   local variables

      integer   ::   unit, io, ierr, logunit

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized. note that data_override_init cannot
!    be called successfully from here (a data_override_mod feature);
!    instead it relies upon a check for previous initialization when
!    subroutine data_override is called.
!---------------------------------------------------------------------
      call fms_init
      call tracer_manager_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=cloudrad_driver_nml, iostat=io)
      ierr = check_nml_error(io,'cloudrad_driver_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=cloudrad_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'cloudrad_driver_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
           write (logunit, nml=cloudrad_driver_nml)

!---------------------------------------------------------------------
!    save the number of shortwave and longwave cloud bands
!---------------------------------------------------------------------
      Cldrad_control%num_sw_cloud_bands = num_sw_cloud_bands
      Cldrad_control%num_lw_cloud_bands = num_lw_cloud_bands

!--------------------------------------------------------------------
!    initialize the modules that are accessed from radiation_driver_mod.
!---------------------------------------------------------------------

      call cloud_spec_init (Exch_ctrl, pref, lonb, latb, axes, Time,   &
                            rad_time_step, Cldrad_control)

      call cloudrad_package_init   (pref, lonb, latb, axes, Time, &
                       Exch_ctrl%donner_meso_is_largescale, Cldrad_control)

!---------------------------------------------------------------------
!    set flag to indicate that module has been successfully initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------

end subroutine cloudrad_driver_init

!######################################################################

subroutine cloudrad_driver (is, ie, js, je,                   &
                            Time, Time_next, Rad_time,        &
                            lat, land, cosz, z_half, z_full,  &
                            r, tsfc, press, pflux, temp, deltaz, &
                            cloudtemp, cloudvapor, clouddeltaz,  &
                            Cldrad_control, &
                            Aerosol, Cld_spec, Model_microphys, &
                            Moist_clouds_block, &
                            crndlw, cmxolw, emrndlw, emmxolw, &
                            camtsw, cldsctsw, cldextsw, cldasymmsw )

!---------------------------------------------------------------------
!    cloudrad_driver obtains cloud radiative properties
!---------------------------------------------------------------------

integer,                 intent(in)         :: is, ie, js, je
type(time_type),              intent(in)  :: Time, Time_next, Rad_time
real, dimension(:,:),         intent(in)             :: lat
real, dimension(:,:),         intent(in)             :: land
real, dimension(:,:,:),       intent(in)             :: z_half, z_full
real, dimension(:,:,:,:),     intent(in)             :: r 
real, dimension(:,:),         intent(in)             :: cosz, tsfc 
real, dimension(:,:,:),       intent(in)             :: press, pflux, &
                                                        temp, deltaz, &
                                                        cloudtemp,  &
                                                        cloudvapor, &
                                                        clouddeltaz
type(cloudrad_control_type),  intent(in)             :: Cldrad_control
type(aerosol_type),           intent(in)             :: Aerosol
type(cld_specification_type), intent(inout)          :: Cld_spec    
type(microphysics_type),      intent(inout)          :: Model_microphys
type(clouds_from_moist_block_type), intent(in)       :: Moist_clouds_block
real, dimension(:,:,:,:),     intent(out)            :: crndlw, camtsw
real, dimension(:,:,:),       intent(out)            :: cmxolw
real, dimension(:,:,:,:,:),   intent(out)            :: emrndlw, emmxolw
real, dimension(:,:,:,:,:),   intent(out)            :: cldsctsw, cldextsw, cldasymmsw

!---------------------------------------------------------------------
!    if the cloud fields are needed, call cloud_spec to retrieve bulk
!    cloud data and place it into a cld_specification_type derived-type 
!    variable Cld_spec and retrieve microphysical data which is returned
!    in microphysics_type variables Lsc_microphys, Meso_microphys and 
!    Cell_microphys and Shallow_microphys, when applicable.
!---------------------------------------------------------------------
type(microphysics_type), allocatable, dimension(:) :: Cloud_microphys

!    Cloud_microphys  microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     cloud schemes, passed through to lower
!                     level routines

integer :: n

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('cloudrad_driver_mod',  &
               'module has not been initialized', FATAL)

!-------------------------------------------------------------------
!    call cloud_spec to retrieve bulk cloud data and place it into
!    a cld_specification_type derived-type variable Cld_spec and
!    retrieve microphysical data which is returned in
!    microphysics_type variables Lsc_microphys, Meso_microphys and 
!    Cell_microphys and Shallow_microphys, when applicable.
!---------------------------------------------------------------------

       call cloud_spec (is, ie, js, je, lat, land, &
                        z_half, z_full, Rad_time,   &
                        press, pflux, temp, cloudtemp, cloudvapor, clouddeltaz, &
                        r, Cldrad_control, Cld_spec, Cloud_microphys, Aerosol, Moist_clouds_block)

!--------------------------------------------------------------------
!    when using the sea-esf radiation, call cloud_radiative_properties
!    to obtain the cloud-radiative properties needed for the radiation 
!    calculation. (these properties are obtained within radiation_calc
!    when executing the original fms radiation code). if these fields 
!    are to be time-averaged, this call is made on all steps; otherwise
!    just on radiation steps.
!--------------------------------------------------------------------

       call cloud_radiative_properties (      &
                     is, ie, js, je, Rad_time, Time, Time_next, &
                     cosz, tsfc, press, pflux, temp, deltaz, &
                     cloudtemp, cloudvapor, clouddeltaz, Cldrad_control, Cld_spec, &
                     Cloud_microphys, Model_microphys, &
                     crndlw, cmxolw, emrndlw, emmxolw, &
                     camtsw, cldsctsw, cldextsw, cldasymmsw)

!---------------------------------------------------------------------
!    deallocate the components of the derived type arrays
!---------------------------------------------------------------------

      do n = 1, size(Cloud_microphys,1)
        call Cloud_microphys(n)%dealloc(Cldrad_control)
      enddo
      deallocate(Cloud_microphys)

!---------------------------------------------------------------------

end subroutine cloudrad_driver

!######################################################################

subroutine cloudrad_driver_end (Cldrad_control)

type(cloudrad_control_type), intent(in) :: Cldrad_control

!----------------------------------------------------------------------
!    cloudrad_driver_end is the destructor
!    for cloudrad_driver_mod.
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('cloudrad_driver_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    wrap up modules specific to the radiation package in use.
!---------------------------------------------------------------------
      call cloudrad_package_end (Cldrad_control)
      call cloud_spec_end (Cldrad_control)

!----------------------------------------------------------------------
!    set initialization status flag.
!----------------------------------------------------------------------
      module_is_initialized = .false.

!----------------------------------------------------------------------

end subroutine cloudrad_driver_end

!######################################################################

             end module cloudrad_driver_mod

