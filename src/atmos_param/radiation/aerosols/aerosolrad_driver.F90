                module aerosolrad_driver_mod

use mpp_mod,               only: input_nml_file
use fms_mod,               only: fms_init, &
                                 mpp_pe, mpp_root_pe, &
                                 open_namelist_file, stdlog, stdout, &
                                 close_file, write_version_number, &
                                 error_mesg, mpp_error, &
                                 check_nml_error, file_exist, &
                                 FATAL, WARNING, NOTE
use time_manager_mod,      only: time_type

use aerosol_types_mod,     only: aerosol_type, &
                                 aerosol_time_vary_type

use aerosol_mod,           only: aerosol_init, aerosol_driver, &
                                 aerosol_time_vary, &
                                 aerosol_endts, &
                                 aerosol_end, &
                                 aerosol_dealloc

use aerosolrad_types_mod,  only: aerosolrad_control_type, &
                                 aerosolrad_diag_type

use aerosolrad_package_mod, only: aerosolrad_package_init,    &
                                  aerosolrad_package_endts, &
                                  aerosolrad_package_time_vary, &
                                  aerosol_radiative_properties, &
                                  aerosolrad_package_end, &
                                  number_of_lw_aerosol_bands

!--------------------------------------------------------------------

implicit none 
private 

!----------------------------------------------------------------------
!------------ version number for this module --------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

!--------------------------------------------------------------------
!------ interfaces -----

public  aerosolrad_driver_init, &
        aerosolrad_driver_time_vary, &
        aerosolrad_driver, &
        aerosolrad_driver_endts, &
        aerosolrad_driver_end, &
        aerosolrad_driver_dealloc

!--------------------------------------------------------------------
!---- namelist ----

logical :: do_swaerosol_forcing = .false.
                                      !  calculating aerosol forcing in
                                      !  shortwave ?
logical :: do_lwaerosol_forcing = .false.
                                      !  calculating aerosol forcing in
                                      !  longwave ?

logical :: override_aerosols_radiation = .false.
                               ! use offline aerosols for radiation 
                               ! calculation
                               ! (via data_override in aerosol_driver)?

namelist /aerosolrad_driver_nml/  do_swaerosol_forcing, &
                                  do_lwaerosol_forcing, &
                                  override_aerosols_radiation

!--------------------------------------------------------------------
!---- private data ----

type(aerosol_time_vary_type) :: Aerosol_rad  ! aerosol_type variable describing
                                             ! the aerosol fields to be seen by
                                             ! the radiation package

logical ::  module_is_initialized = .false. ! module initialized?

!---------------------------------------------------------------------

                         contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!######################################################################

subroutine aerosolrad_driver_init (lonb, latb, kmax,   &
                                   num_sw_bands, num_lw_bands, &
                                   Aerosolrad_control, &
                                   aerosol_names, aerosol_family_names)

real, dimension(:,:),       intent(in)    :: lonb, latb
integer,                    intent(in)    :: kmax
integer,                    intent(in)    :: num_sw_bands, &
                                             num_lw_bands
type(aerosolrad_control_type), intent(inout) :: Aerosolrad_control
character(len=64), pointer, intent(out)   :: aerosol_names(:), &
                                             aerosol_family_names(:)

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
      integer  ::   unit, io, ierr, logunit

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module have been initialized
!---------------------------------------------------------------------
      call fms_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=aerosolrad_driver_nml, iostat=io)
      ierr = check_nml_error(io,'aerosolrad_driver_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosolrad_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'aerosolrad_driver_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
           write (logunit, nml=aerosolrad_driver_nml)

!---------------------------------------------------------------------
!    define control variables indicating whether the aerosol forcings
!    should be calculated, also set the size of the output for
!    aerosol forcing
!---------------------------------------------------------------------
      Aerosolrad_control%do_lwaerosol_forcing = do_lwaerosol_forcing
      Aerosolrad_control%do_swaerosol_forcing = do_swaerosol_forcing
      Aerosolrad_control%indx_lwaf = 1
      Aerosolrad_control%indx_swaf = 1
      if (do_lwaerosol_forcing) Aerosolrad_control%indx_lwaf = Aerosolrad_control%indx_lwaf + 1
      if (do_swaerosol_forcing) Aerosolrad_control%indx_swaf = Aerosolrad_control%indx_swaf + 1

!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize aerosol modules    
!-----------------------------------------------------------------------
      call aerosol_init (lonb, latb, Aerosol_rad, &
                         aerosol_names, aerosol_family_names)

      call aerosolrad_package_init (kmax, aerosol_names, lonb, latb, &
                                    Aerosolrad_control)

!---------------------------------------------------------------------
!    save the number of shortwave and longwave aerosols
!    NOTE: number of LW aerosols is determined by aerosolrad package
!---------------------------------------------------------------------
      Aerosolrad_control%num_sw_aerosol_bands = num_sw_bands
      call number_of_lw_aerosol_bands (Aerosolrad_control%num_lw_aerosol_bands)

!---------------------------------------------------------------------

      module_is_initialized = .true.

!---------------------------------------------------------------------

end subroutine aerosolrad_driver_init

!######################################################################

subroutine aerosolrad_driver_time_vary (Rad_time, Aerosolrad_control)
type(time_type),               intent(in) :: Rad_time
type(aerosolrad_control_type), intent(in) :: Aerosolrad_control

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('aerosolrad_driver_mod',  &
               'module has not been initialized', FATAL)

!-------------------------------------------------------------------
      call aerosol_time_vary (Rad_time, Aerosol_rad)
      call aerosolrad_package_time_vary (Rad_time, Aerosolrad_control)
!-------------------------------------------------------------------

end subroutine aerosolrad_driver_time_vary

!######################################################################

subroutine aerosolrad_driver_endts

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('aerosolrad_driver_mod',  &
               'module has not been initialized', FATAL)

!-------------------------------------------------------------------
      call aerosol_endts (Aerosol_rad)
      call aerosolrad_package_endts
!-------------------------------------------------------------------

end subroutine aerosolrad_driver_endts

!######################################################################

subroutine aerosolrad_driver (is, ie, js, je, Rad_time, &
                              fracday, phalf, pflux, &
                              aerosolrelhum, deltaz, r, &
                              do_cmip_sw_diagnostics, Aerosolrad_control, &
                              Aerosolrad_diags, Aerosol, &
                              extinction, aerooptdep, aerooptdep_volc, &
                              aeroasymfac, aerosctopdep, aeroextopdep)
                                     
!-------------------------------------------------------------------
integer,                       intent(in)  :: is, ie, js, je
type(time_type),               intent(in)  :: Rad_time
real, dimension(:,:),          intent(in)  :: fracday
real, dimension(:,:,:),        intent(in)  :: phalf, pflux, &
                                              aerosolrelhum, deltaz
real, dimension(:,:,:,:),      intent(in)  :: r
logical,                       intent(in)  :: do_cmip_sw_diagnostics
type(aerosolrad_control_type), intent(in)  :: Aerosolrad_control
type(aerosolrad_diag_type),    intent(out) :: Aerosolrad_diags
type(aerosol_type),            intent(out) :: Aerosol
real, dimension(:,:,:),        intent(out) :: extinction
real, dimension(:,:,:,:),      intent(out) :: aerooptdep, aerooptdep_volc, &
                                              aeroasymfac, aerosctopdep, aeroextopdep
!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('aerosolrad_driver_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    if the aerosol fields are needed as input to the radiation package,
!    call aerosol_driver to access the aerosol data and place it into 
!    an aerosol_type derived-type variable Aerosol.
!---------------------------------------------------------------------
      call aerosol_driver (is, js, Rad_time, r, phalf, pflux, &
                           Aerosol_rad, Aerosol, override_aerosols_radiation)

!---------------------------------------------------------------------
!    obtain the aerosol optical properties for longwave and shortwave
!---------------------------------------------------------------------
      call aerosol_radiative_properties (is, ie, js, je, &
                                         Rad_time, fracday, &
                                         pflux, aerosolrelhum, deltaz, &
                                         do_cmip_sw_diagnostics, &
                                         Aerosolrad_control, &
                                         Aerosolrad_diags, Aerosol, &
                                         extinction, aerooptdep, aerooptdep_volc, &
                                         aeroasymfac, aerosctopdep, aeroextopdep)

!---------------------------------------------------------------------

end subroutine aerosolrad_driver

!######################################################################

subroutine aerosolrad_driver_end (Aerosolrad_control)

type(aerosolrad_control_type), intent(in) :: Aerosolrad_control

!---------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('aerosolrad_driver_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------
      call aerosolrad_package_end (Aerosolrad_control)
      call aerosol_end (Aerosol_rad)
!---------------------------------------------------------------------

      module_is_initialized = .false.

!---------------------------------------------------------------------

end subroutine aerosolrad_driver_end

!######################################################################

subroutine aerosolrad_driver_dealloc (Aerosol, Aerosolrad_diags)

type(aerosol_type),         intent(inout) :: Aerosol
type(aerosolrad_diag_type), intent(inout) :: Aerosolrad_diags

!---------------------------------------------------------------------
! deallocate the local instanance of several aerosol types
! the data structures were allocated by aerosolrad_driver
!---------------------------------------------------------------------

      call aerosol_dealloc (Aerosol)
     !call aerosolrad_diag_dealloc (Aerosolrad_diags)
      call Aerosolrad_diags%dealloc

!-------------------------------------------------------------------

end subroutine aerosolrad_driver_dealloc

!######################################################################

                end module aerosolrad_driver_mod

