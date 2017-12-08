                     module longwave_driver_mod
!
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Code to set up longwave radiation calculation
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
! 

use mpp_mod,            only: input_nml_file
use fms_mod,            only: open_namelist_file, fms_init, &
                              mpp_pe, mpp_root_pe, stdlog, &
                              file_exist, write_version_number, &
                              check_nml_error, error_mesg, &
                              FATAL, close_file

! shared radiation package modules:

use radiation_driver_types_mod, only: radiation_control_type

use aerosolrad_types_mod, only: aerosolrad_control_type

!  longwave radiation package modules:

use sealw99_mod,        only: sealw99_init, sealw99_time_vary, sealw99, &
                              sealw99_endts, sealw99_end, &
                              lw_table_type, &
                              longwave_number_of_bands => sealw99_number_of_bands, &
                              longwave_get_tables => sealw99_get_tables

use longwave_types_mod, only: lw_output_type, lw_diagnostics_type, &
                              assignment(=)

use radiative_gases_mod,       only: get_longwave_gas_flag
use radiative_gases_types_mod, only: radiative_gases_type

!------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!    longwave_driver_mod is the driver for the longwave radiation
!    component of the sea_esf_rad radiation package.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------

public      &
   longwave_driver_init, longwave_driver_time_vary, longwave_driver,   &
   longwave_driver_endts, longwave_driver_end

! inherited from sealw99_mod
public  longwave_number_of_bands, longwave_get_tables, &
        lw_table_type

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16) :: lwform= 'sealw99'
logical :: lw_cs = .false.        !LGS: .True. for longwave-cloud-radiative effects turned off. 

integer :: clear_min_lw = 0    !mtp: If lw_cs .EQ. .true., clear_min_lw specifies the minimum layer at which
                               !     clouds are turned off. (only works if clear_max_lw .ne. 0 as well)
integer :: clear_max_lw = 0    !mtp: If lw_cs .EQ. .true., clear_max_lw specifies the maximum layer at which
                               !     clouds are turned off. (only works if clear_min_lw .ne. 0 as well)

logical :: no_lw_cloud_heating_at_levs = .false. ! mtp: Turns of shortwave cloud heating 
                                                 ! from level no_cloud_min to level
                                                 ! no_cloud_max if .true.

integer :: no_heating_min_lw = 1           ! level=1 for top layer         
integer :: no_heating_max_lw = 1           ! level= 32/48 is bottom layer (in 32/48 layer model) 
 



namelist / longwave_driver_nml /   lwform, lw_cs,                        &
                                   no_lw_cloud_heating_at_levs,          &
                                   no_heating_min_lw, no_heating_max_lw, &
                                   clear_min_lw, clear_max_lw                

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

logical :: module_is_initialized =  .false.   ! module initialized ?
logical :: do_sealw99 = .false.               ! sealw99 parameter-
                                              ! ization active ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!                    PUBLIC SUBROUTINES
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! <SUBROUTINE NAME="longwave_driver_init">
!  <OVERVIEW>
!   longwave_driver_init is the constructor for longwave_driver_mod
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine initializes longwave radiation package
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_driver_init (pref)
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!   array containing two reference pressure profiles [pascals]
!  </IN>
! </SUBROUTINE>
!
subroutine longwave_driver_init (pref)
 
!---------------------------------------------------------------------
!    longwave_driver_init is the constructor for longwave_driver_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:,:),         intent(in) :: pref

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 [ Pa ]
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!  local variables

      integer     :: unit, ierr, io, logunit

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=longwave_driver_nml, iostat=io)
      ierr = check_nml_error(io,'longwave_driver_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=longwave_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'longwave_driver_nml')
        end do
10      call close_file (unit)
      endif
#endif
 
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=longwave_driver_nml)

!--------------------------------------------------------------------
!    determine if valid specification of lw radiation has been made.
!    if optional packages are provided at some later time, this is where
!    the choice of package will be made.
!---------------------------------------------------------------------
      if (trim(lwform) == 'sealw99') then
        do_sealw99 = .true.
        call sealw99_init ( pref, &
                            get_longwave_gas_flag('h2o'), &
                            get_longwave_gas_flag('o3'),  &
                            get_longwave_gas_flag('ch4'), &
                            get_longwave_gas_flag('n2o'), &
                            get_longwave_gas_flag('co2'), &
                            get_longwave_gas_flag('co2_10um'), &
                            get_longwave_gas_flag('cfc') )
      else
        call error_mesg ( 'longwave_driver_mod', &
                 'invalid longwave radiation form specified', FATAL)
      endif

!---------------------------------------------------------------------
!    set flag indicating successful initialization of module.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------


end  subroutine longwave_driver_init


!#####################################################################

subroutine longwave_driver_time_vary (Rad_gases_tv)
 
!-------------------------------------------------------------------- 
type(radiative_gases_type), intent(inout) :: Rad_gases_tv

      call sealw99_time_vary (Rad_gases_tv%use_ch4_for_tf_calc, &
                              Rad_gases_tv%use_n2o_for_tf_calc, &
                              Rad_gases_tv%use_co2_for_tf_calc, &
                              Rad_gases_tv%ch4_for_tf_calc, &
                              Rad_gases_tv%n2o_for_tf_calc, &
                              Rad_gases_tv%co2_for_tf_calc)
 
end subroutine longwave_driver_time_vary      

        
!#####################################################################

subroutine longwave_driver_endts (Rad_gases_tv)
          
type(radiative_gases_type), intent(in) :: Rad_gases_tv
 

     call sealw99_endts   ! does nothing


end subroutine longwave_driver_endts 



!#####################################################################
! <SUBROUTINE NAME="longwave_driver">
!  <OVERVIEW>
!   Subroutine to set up and execute longwave radiation calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   longwave_driver allocates and initializes longwave radiation out-
!    put variables and selects an available longwave radiation param-
!    eterization, executes it, and then returns the output fields to 
!    sea_esf_rad_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_driver (Atmos_input, Rad_gases, &
!                         emrndlw, emmxolw, crndlw, cmxolw, &
!                         aerooptdep, aerooptdep_volc, &
!                         Lw_output, Lw_diagnostics)
!
!  </TEMPLATE>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!  </IN>
!  <IN NAME="emrndlw" TYPE="real">
!   cloud emissivity for random overlap clouds
!   by longwave band and profile
!  </IN>
!  <IN NAME="emmxolw" TYPE="real">
!   cloud emissivity for maximum overlap clouds
!   by longwave band and profile
!  </IN>
!  <IN NAME="crndlw" TYPE="real">
!   cloud amount for random overlap clouds
!  </IN>
!  <IN NAME="cmxolw" TYPE="real">
!   cloud amount for maximum overlap clouds
!  </IN>
!  <IN NAME="aerooptdep" TYPE="real">
!   Aerosol optical depth in model layers
!  </IN>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                   radiation output data
!  </INOUT>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!  </INOUT>
! </SUBROUTINE>
!
subroutine longwave_driver (press, pflux, temp, tflux, rh2o, deltaz, &
                            Rad_gases, emrndlw, emmxolw, crndlw, cmxolw, &
                            aerooptdep, aerooptdep_volc, &
                            flag_stoch, Rad_control, &
                            do_lwaerosol, volcanic_lw_aerosols, &
                            Lw_output, Lw_diagnostics)

!--------------------------------------------------------------------
!    longwave_driver allocates and initializes longwave radiation out-
!    put variables and selects an available longwave radiation param-
!    eterization, executes it, and then returns the output fields to 
!    sea_esf_rad_mod.
!--------------------------------------------------------------------

real, dimension(:,:,:),       intent(in)     :: press, pflux, temp, &
                                                tflux, rh2o, deltaz
type(radiative_gases_type),   intent(inout)  :: Rad_gases   
real, dimension(:,:,:,:,:),   intent(in)     :: emrndlw, emmxolw
real, dimension(:,:,:,:),     intent(in)     :: crndlw
real, dimension(:,:,:),       intent(in)     :: cmxolw
real, dimension(:,:,:,:),     intent(in)     :: aerooptdep, aerooptdep_volc
integer,                      intent(in)     :: flag_stoch
type(radiation_control_type),  intent(in)    :: Rad_control
logical,                       intent(in)    :: do_lwaerosol
logical,                       intent(in)    :: volcanic_lw_aerosols
type(lw_output_type), dimension(:),  intent(inout)  :: Lw_output
type(lw_diagnostics_type),    intent(inout)  :: Lw_diagnostics

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      press          pressure at model levels [pascals]
!      pflux          pressure at model layer interfaces [pascals]
!      temp           mean temperature of model layer [kelvin]
!      tflux          temperature at model layer interfaces [kelvin]
!      rh2o           mixing ratio for water vapor
!      deltaz         model layer thickness in meters
!      Rad_gases      radiative_gases_type variable containing the radi-
!                     ative gas input fields needed by the radiation 
!                     package
!      emrndlw        cloud emissivity for random overlap clouds by band
!      emmxolw        cloud emissivity for maximum overlap clouds by band
!      crndlw         cloud amount for random overlap clouds
!      cmxolw         cloud amount for maximum overlap clouds
!      aerooptdep     aerosol optical depth in model layers
!
!   intent(inout) variables:
!
!      Lw_output      lw_output_type variable containing longwave 
!                     radiation output data 
!      Lw_diagnostics lw_diagnostics_type variable containing diagnostic
!                     longwave output used by the radiation diagnostics
!                     module
!  
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables

      logical :: do_lwaerosol_forcing
      logical :: calc_includes_aerosols
      integer :: ix, jx, kx  ! dimensions of current physics window
      integer :: indx

      real    :: crndlw_local(size(crndlw,1),size(crndlw,2),size(crndlw,3),size(crndlw,4))
      real    :: cmxolw_local(size(cmxolw,1),size(cmxolw,2),size(cmxolw,3))

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_driver_mod',   &
          'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    call longwave_driver_alloc to allocate component arrays of a
!    lw_output_type variable.
!----------------------------------------------------------------------
      ix = size(press,1)
      jx = size(press,2)
      kx = size(press,3)

      do indx = 1, size(Lw_output,1)
         call Lw_output(indx)%alloc (ix, jx, kx, Rad_control%do_totcld_forcing)
      enddo

      if (size(Lw_output,1) .gt. 1) then
         do_lwaerosol_forcing = .true.
      else
         do_lwaerosol_forcing = .false.
      endif

!--------------------------------------------------------------------
!    calculate the longwave radiative heating rates and fluxes.
!--------------------------------------------------------------------
      if (do_sealw99) then

          !LGS
          crndlw_local=crndlw
          cmxolw_local=cmxolw
          if (lw_cs) then
             if (clear_min_lw .EQ. 0 .OR. clear_max_lw .EQ. 0) then
               crndlw_local=0.0
               cmxolw_local=0.0
             else
               crndlw_local(:,:,clear_min_lw:clear_max_lw,:)=0.0
               cmxolw_local(:,:,clear_min_lw:clear_max_lw)=0.0
             end if
          endif
             
!--------------------------------------------------------------------
!    call sealw99 to use the simplified-exchange-approximation (sea)
!    parameterization.
!----------------------------------------------------------------------
         if (do_lwaerosol_forcing) then
           if (do_lwaerosol) then
             calc_includes_aerosols = .false.
           else
             calc_includes_aerosols = .true.
           endif

!----------------------------------------------------------------------
!    call sealw99 with aerosols (if model is being run without) and 
!    without aerosols (if model is being run with). save the radiation
!    fluxes to Lw_output_ad (which does not feed back into the model),
!    but which may be used to define the aerosol forcing.
!----------------------------------------------------------------------
           indx = size(Lw_output,1)
           call sealw99 ( &
                     press, pflux, temp, tflux, rh2o, deltaz, &
                     Rad_gases%qo3, Rad_gases%rrvco2,      &
                     Rad_gases%rrvf11, Rad_gases%rrvf12,   &
                     Rad_gases%rrvf113, Rad_gases%rrvf22,  &
                     emrndlw, emmxolw, crndlw_local, cmxolw_local, &
                     aerooptdep, aerooptdep_volc, Lw_output(indx), &
                     Lw_diagnostics, flag_stoch, &
                     Rad_control%do_totcld_forcing, calc_includes_aerosols, &
                     volcanic_lw_aerosols)
!-------------------------------------------------------------------------------
 ! mtp: replacing lw-all-sky heating rates by clear sky heating rates from
 ! level no_cloud_min to level no_cloud_max
 ! level 1=top layer, level 32 = surface layer (for 32 layer atmosphere)
 ! levl 20 = 663 hPa

           if (no_lw_cloud_heating_at_levs) then
             Lw_output(indx)%heatra(:,:,no_heating_min_lw:no_heating_max_lw) = &
              Lw_output(indx)%heatracf(:,:,no_heating_min_lw:no_heating_max_lw)
           end if
!--------------------------------------------------------------------------------
         endif
 
!----------------------------------------------------------------------
!    standard call, where radiation output feeds back into the model.
!----------------------------------------------------------------------
           call sealw99 ( &
                     press, pflux, temp, tflux, rh2o, deltaz, &
                     Rad_gases%qo3, Rad_gases%rrvco2,      &   
                     Rad_gases%rrvf11, Rad_gases%rrvf12,   &   
                     Rad_gases%rrvf113, Rad_gases%rrvf22,  &
                     emrndlw, emmxolw, crndlw_local, cmxolw_local, &
                     aerooptdep, aerooptdep_volc, Lw_output(1),  &
                     Lw_diagnostics, flag_stoch, &
                     Rad_control%do_totcld_forcing, do_lwaerosol, &
                     volcanic_lw_aerosols)
 
!-------------------------------------------------------------------------------
 ! mtp: replacing lw-all-sky heating rates by clear sky heating rates from
 ! level no_cloud_min to level no_cloud_max
 ! level 1=top layer, level 32 = surface layer (for 32 layer atmosphere)
 ! levl 20 = 663 hPa

if (no_lw_cloud_heating_at_levs) then
  Lw_output(1)%heatra(:,:,no_heating_min_lw:no_heating_max_lw) = Lw_output(1)%heatracf(:,:,no_heating_min_lw:no_heating_max_lw)
end if
!--------------------------------------------------------------------------------
      else

!--------------------------------------------------------------------
!    at the current time sealw99 is the only longwave parameterization 
!    available.
!----------------------------------------------------------------------
        call error_mesg ('longwave_driver_mod', &
         'invalid longwave radiation parameterization selected', FATAL)
      endif

!---------------------------------------------------------------------

end subroutine longwave_driver

!#####################################################################
! <SUBROUTINE NAME="longwave_driver_end">
!  <OVERVIEW>
!   Subroutine to end longwave calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine end longwave calculation
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_driver_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine longwave_driver_end                  

!--------------------------------------------------------------------
!    longwave_driver_end is the destructor for longwave_driver_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_driver_mod',   &
          'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    call sealw99_end to close sealw99_mod.
!-------------------------------------------------------------------
      if (do_sealw99) then
        call sealw99_end
      endif

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

end subroutine longwave_driver_end                  

!###################################################################

                end module longwave_driver_mod

