                     module shortwave_driver_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
!
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
! 
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
!
! <OVERVIEW>
!  Code to carry out shortwave calculation.
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes, prepares, and ends shortwave radiation calculation.
!  This code is called by sea_esf_rad.f90 and calls shortwave subroutines
!  to do shortwave flux calculation.
! </DESCRIPTION>
!

!   shared modules:

use mpp_mod,              only: input_nml_file
use fms_mod,              only: open_namelist_file, fms_init, &
                                mpp_pe, mpp_root_pe, stdlog, &
                                file_exist, write_version_number, &
                                check_nml_error, error_mesg, &
                                FATAL, NOTE, close_file
use time_manager_mod,     only: time_manager_init, time_type, &
                                set_date, get_date, print_date, &
                                assignment(=), operator(-), &
                                operator(>), operator(+)
use diag_manager_mod,     only: diag_manager_init, get_base_time


!   shared radiation package modules:
 
use radiation_driver_types_mod, only: radiation_control_type, &
                                      astronomy_type

use esfsw_driver_mod,     only: esfsw_driver_init, swresf,   &
                                esfsw_driver_end, &
                                shortwave_number_of_bands => esfsw_number_of_bands, &
                                esfsw_solar_flux

!  radiation package modules:

use shortwave_types_mod,  only: sw_output_type, assignment(=)

use radiative_gases_types_mod, only: radiative_gases_type

use solar_data_driver_mod, only: solar_data_driver_init, &
                                 solar_data_driver_time_vary, &
                                 solar_data_driver_end

!-------------------------------------------------------------------

implicit none
private

!------------------------------------------------------------------
!    shortwave_driver_mod is the driver for shortwave radiation 
!    component of the sea_esf_rad radiation package.
!-----------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module  -------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public        &
          shortwave_driver_init , shortwave_driver,    & 
          shortwave_driver_end, &
          shortwave_driver_time_vary, &
          get_solar_flux_by_band, &
          get_solar_constant

public    shortwave_number_of_bands


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)   :: swform = '    '
logical             :: do_cmip_diagnostics = .false.
!logical             :: calculate_volcanic_sw_heating = .false.

real    :: solar_constant = 1365.0    !  annual mean solar flux at top 
                                      !  of atmosphere [ W/(m**2) ]

logical :: time_varying_solar_constant = .false.
                                      !  solar_constant is to vary with
                                      !  time ?

integer, dimension(6) :: solar_dataset_entry = (/ 1, 1, 1, 0, 0, 0 /)
                                      ! time in solar data set corresp-
                                      ! onding to model initial time
                                      ! (yr, mo, dy, hr, mn, sc)

logical :: sw_cs = .false.            ! Turns of shortwave cloud radiative effects if .true.

integer :: clear_min_sw = 0    !mtp: If sw_cs .EQ. .true., clear_min_lw specifies the minimum layer at which
                               !     clouds are turned off. (only works if clear_max_lw .ne. 0 as well)
integer :: clear_max_sw = 0    !mtp: If sw_cs .EQ. .true., clear_max_lw specifies the maximum layer at which
                               !     clouds are turned off. (only works if clear_min_lw .ne. 0 as well)



logical :: no_sw_cloud_heating_at_levs = .false. ! mtp: Turns of shortwave cloud heating 
                                                 ! from level no_cloud_min to level
                                                 ! no_cloud_max if .true.

integer :: no_heating_min_sw = 1        ! level=1 for top layer         
integer :: no_heating_max_sw = 1        ! level=32/48 for bottom layer (in 32/48 layer model) 


 
namelist / shortwave_driver_nml /    do_cmip_diagnostics, &
!                                    calculate_volcanic_sw_heating, &
                                     swform, &
                                     time_varying_solar_constant, &
                                     solar_dataset_entry, &
                                     solar_constant, &
                                     sw_cs, &
                                     no_sw_cloud_heating_at_levs, &
                                     no_heating_min_sw, &
                                     no_heating_max_sw, &
                                     clear_min_sw, &
                                     clear_max_sw  
                                     

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

! solar data needed by the esfsw code
real, allocatable :: solflxband(:)
real              :: solar_constant_used
logical           :: solflxband_initialized
! reference bands may be needed by chemistry
real, allocatable :: solflxbandref(:)

logical         :: negative_offset
type(time_type) :: Solar_offset

logical :: module_is_initialized = .false.  ! module initialized ?


!-------------------------------------------------------------------
!-------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! <SUBROUTINE NAME="shortwave_driver_init">
!  <OVERVIEW>
!   Code that initializes shortwave radiation calculation.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Initialize utilities and radiation utilities if necessary. They
!   should have been initialized in the radiation initialiation subroutine
!   in the sea_esf_rad.f90. The code then reads in input.nml namelist
!   and logs input parameters to logfile. It uses lhsw or esfsw package
!   depends on namelist parameter. Initializes apropriate shortwave
!   package subroutines and set up the initialize parameter.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call shortwave_driver_init (latb, pref)
!  </TEMPLATE>
!  <IN NAME="latb" TYPE="real">
!   2d array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="pref" TYPE="real">
!   An array containing two reference pressure profiles [pascals]
!  </IN>
! </SUBROUTINE>
subroutine shortwave_driver_init (Rad_control)

!---------------------------------------------------------------------
!    shortwave_driver_init is the constructor for shortwave_driver_mod.
!---------------------------------------------------------------------
type(radiation_control_type), intent(inout) :: Rad_control


!---------------------------------------------------------------------
!  local variables:

      integer   :: unit, io, ierr, logunit
      integer   :: nbands

      type(time_type) :: Time_init, Solar_entry
!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!                                
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!-------------------------------------------------------------------
      if (module_is_initialized) return

!-------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!-------------------------------------------------------------------
      call fms_init
      call diag_manager_init
      call time_manager_init

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=shortwave_driver_nml, iostat=io)
      ierr = check_nml_error(io,"shortwave_driver_nml")
#else
!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=shortwave_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'shortwave_driver_nml')
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
                       write (logunit, nml=shortwave_driver_nml)

!--------------------------------------------------------------------
!    define logicals specifying the sw package in use as an element
!    of a shortwave_control_type variable usable by other radiation-
!    related modules. initialize the modules associated with the chosen
!    sw package.
!---------------------------------------------------------------------
      if (trim(swform) == 'esfsw99') then
        call esfsw_driver_init
        call shortwave_number_of_bands(nbands)
        ! save reference solar data
        allocate (solflxbandref(nbands))
        call esfsw_solar_flux (solflxbandref, ref=.true.)
      else
        call error_mesg ( 'shortwave_driver_mod',   &
        'improper specification of desired shortwave parameterization',&
                                                               FATAL)
      endif

!---------------------------------------------------------------------
!    save the logical indicating the need to generate cmip aerosol
!    diagnostics and mark it as initialized.
!---------------------------------------------------------------------
      Rad_control%do_cmip_sw_diagnostics = do_cmip_diagnostics

!---------------------------------------------------------------------
!    set up solar data
!---------------------------------------------------------------------

      call solar_data_driver_init (nbands, ierr)

      ! if the error code returned is non-zero then
      ! the data file did not exist
      if (ierr /= 0) then
         if (time_varying_solar_constant) then
            call error_mesg ('shortwave_driver_mod', &
            'desired solar_spectral_data input file is not present', FATAL)
         endif
      endif
      
      allocate (solflxband(nbands))
      solflxband_initialized = .false.

      if (time_varying_solar_constant) then

         Time_init = get_base_time()  ! this is a bug - should be Time_init

         call compute_time_offset ( Time_init, solar_dataset_entry, &
                                    Solar_entry, Solar_offset, negative_offset )

         call error_mesg ('shortwave_driver_mod', &
               'Solar data is varying in time', NOTE)
         call print_date (Solar_entry, &
               str='Data from solar timeseries at time: ')
      else

         if (solar_dataset_entry(1) == 1 .and. &
             solar_dataset_entry(2) == 1 .and. &
             solar_dataset_entry(3) == 1 .and. &
             solar_dataset_entry(4) == 0 .and. &
             solar_dataset_entry(5) == 0 .and. &
             solar_dataset_entry(6) == 0 ) then
               ! no solar date
               ! use solar constant from namelist and reference solar by band
              !call esfsw_solar_flux_init (solar_constant, solflxbandref)
               solar_constant_used = solar_constant
               solflxband = solflxbandref
               solflxband_initialized = .true.
               call error_mesg ('shortwave_driver_mod', &
                    'Solar data is fixed in time at nml value', NOTE)
         else
               ! use nml solar date
               ! convert date to time_type variable
               Solar_entry = set_date (solar_dataset_entry(1), &
                                       solar_dataset_entry(2), &
                                       solar_dataset_entry(3), &
                                       solar_dataset_entry(4), &
                                       solar_dataset_entry(5), &
                                       solar_dataset_entry(6))
               call error_mesg ('shortwave_driver_mod', &
                            'Solar data is fixed in time', NOTE)
               call print_date (Solar_entry ,    &
                     str='Data used in this experiment is from solar &
                          &timeseries at time:')
               ! define time to be used for solar input data
               call solar_data_driver_time_vary (Solar_entry, solar_constant_used, solflxband)
               solflxband_initialized = .true.
              !call esfsw_solar_flux_init (solar_constant_used, solflxband)
         endif
      endif


!-------------------------------------------------------------------
!    set flag indicating successful initialization of module.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------


end subroutine shortwave_driver_init



!###########################################################

subroutine shortwave_driver_time_vary (Rad_time)
type(time_type),    intent(in)  :: Rad_time

type(time_type) :: Solar_time

      if (time_varying_solar_constant) then

         ! compute offset from radiation time
         ! may remove this option in the future
         if (negative_offset) then
            Solar_time = Rad_time - Solar_offset
         else
            Solar_time = Rad_time + Solar_offset
         endif

        call solar_data_driver_time_vary (Solar_time, solar_constant_used, solflxband)
       !call esfsw_solar_flux_init (solar_constant_used, solflxband)
        solflxband_initialized = .true.
      endif

end subroutine shortwave_driver_time_vary

!###########################################################
! <SUBROUTINE NAME="shortwave_driver">
!  <OVERVIEW>
!   Code that deploys shortwave radiation calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!    shortwave_driver initializes shortwave radiation output variables, 
!    determines if shortwave radiation is present in the current physics
!    window, selects one of the available shortwave parameterizations,
!    executes it, and returns the output fields to sea_esf_rad_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call shortwave_driver (Atmos_input, Surface, Astro, &
!                           aeroasymfac, aerosctopdep, aeroextopdep,     &
!                           Rad_gases, camtsw, cldsct, cldext, cldasymm, &
!                           Sw_output)
!  </TEMPLATE>
!   <IN NAME="Astro" TYPE="astronomy_type">
!     Astronomy_type variable containing the astronomical
!     input fields on the radiation grid  
!   </IN>
!   <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!     Radiative_gases_type variable containing the radiative 
!     gas input fields on the radiation grid 
!   </IN>
!   <IN NAME="camtsw" TYPE="real">
!    Cloud amount for shortwave clouds. If stochastic clouds is implemented
!    then cloud amount by shortwave band.
!   </IN>
!   <IN NAME="cldext" TYPE="real">
!    Cloud extinction parameter
!   </IN>
!   <IN NAME="cldsct" TYPE="real">
!    Cloud single scattering albedo
!   </IN>
!   <IN NAME="cldasymm" TYPE="real">
!    Cloud asymmetric parameter
!   </IN>
!   <INOUT NAME="Sw_output" TYPE="sw_output_type">
!     The shortwave radiation calculation result
!   </INOUT>
! </SUBROUTINE>
subroutine shortwave_driver (press, pflux, temp, rh2o, &
                             deltaz, asfc_vis_dir, asfc_nir_dir, &
                             asfc_vis_dif, asfc_nir_dif, Astro,   &
                             aeroasymfac, aerosctopdep, aeroextopdep,       &
                             Rad_gases, camtsw, cldsct, cldext, cldasymm,   &
                             flag_stoch, Rad_control, do_swaerosol, Sw_output )

!---------------------------------------------------------------------
!    shortwave_driver initializes shortwave radiation output variables, 
!    determines if shortwave radiation is present in the current physics
!    window, selects one of the available shortwave parameterizations,
!    executes it, and returns the output fields to sea_esf_rad_mod.
!---------------------------------------------------------------------

real, dimension(:,:,:),          intent(in)    :: press, pflux, temp, rh2o, deltaz
real, dimension(:,:),            intent(in)    :: asfc_vis_dir, asfc_nir_dir, &
                                                  asfc_vis_dif, asfc_nir_dif
type(astronomy_type),            intent(in)    :: Astro           
type(radiative_gases_type),      intent(in)    :: Rad_gases   
real,dimension(:,:,:,:),         intent(in)    :: aeroasymfac, aerosctopdep, aeroextopdep
real, dimension(:,:,:,:),        intent(in)    :: camtsw
real, dimension(:,:,:,:,:),      intent(in)    :: cldsct, cldext, cldasymm
integer,                         intent(in)    :: flag_stoch
type(radiation_control_type),    intent(in)    :: Rad_control
logical,                         intent(in)    :: do_swaerosol
type(sw_output_type), dimension(:), intent(inout) :: Sw_output

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      press          pressure at model levels [pascals]
!      pflux          pressure at model layer interfaces [pascals]
!      temp           mean temperature of model layer [kelvin]
!      rh2o           mixing ratio for water vapor
!      deltaz         model layer thickness in meters
!      Astro          astronomy_type variable containing the astronom-
!                     ical input fields needed by the radiation package
!      Rad_gases      radiative_gases_type variable containing the radi-
!                     ative gas input fields needed by the radiation 
!                     package
!      aeroasymfac, aerosctopdep, aeroextopdep
!                     aerosol radiative properties
!      camtsw         cloud amount for shortwave clouds,
!                     if stochastic clouds is implemented then
!                     cloud amount by shortwave band
!      cldsct, cldext, cldasymm
!                     cloud radiative properties
!
!   intent(out) variables:
!
!      Sw_output      sw_output_type variable containing shortwave 
!                     radiation output data 
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      logical  :: skipswrad
      logical  :: use_aero
      logical  :: do_swaerosol_forcing
      integer  :: naerosol_optical
      integer  :: i, j, indx
      integer  :: ix, jx, kx
      real     :: aerozero(size(aeroasymfac,1), &
                           size(aeroasymfac,2), &
                           size(aeroasymfac,3), &
                           size(aeroasymfac,4))

      real     :: camtsw_local(size(camtsw,1),size(camtsw,2),size(camtsw,3),size(camtsw,4))
!---------------------------------------------------------------------
!   local variables:
!
!      skipswrad    bypass calling sw package because sun is not 
!                   shining any where in current physics window ?
!      do_swaerosol_forcing  if the output (type) array contains
!                   more than one-dimension then aerosol forcing
!                   is to be computed
!      ix,jx,kx     dimensions of current physics window
!      i,j          do-loop indices
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('shortwave_driver_mod',   &
             'module has not been initialized', FATAL )
      endif

!----------------------------------------------------------------------
!    call shortwave_output_alloc to initialize shortwave fluxes and 
!    heating rates.
!--------------------------------------------------------------------
      ix = size(press,1)
      jx = size(press,2)
      kx = size(press,3)

      do indx = 1, size(Sw_output,1)
         call Sw_output(indx)%alloc (ix, jx, kx, Rad_control%do_totcld_forcing)
      enddo

      if (size(Sw_output,1) .gt. 1) then
         do_swaerosol_forcing = .true.
      else
         do_swaerosol_forcing = .false.
      endif

!--------------------------------------------------------------------
!    determine when the no-sun case exists at all points within the 
!    physics window and bypass the sw radiation calculations for that 
!    window. for do_annual_mean or do_daily_mean, only one cosz in a
!    model row need be tested, since all points in i have the same 
!    zenith angle.
!--------------------------------------------------------------------
      skipswrad = .true.
      do j=1,jx        
        if ( Astro%cosz(1,j) > 0.0 ) skipswrad = .false.
        if (Rad_control%do_diurnal) then
          do i = 2,ix         
            if (Astro%cosz(i,j) > 0.0 )  then
              skipswrad = .false.
              exit
            endif
          end do
        endif
      end do

!--------------------------------------------------------------------
!    for aerosol optical depth diagnostics, swresf must be called
!    on all radiation steps.
!--------------------------------------------------------------------
      if (do_cmip_diagnostics)  skipswrad = .false.

!--------------------------------------------------------------------
!    if the sun is shining nowhere in the physics window allocate
!    output fields which will be needed later, set them to a flag
!    value and return.
!--------------------------------------------------------------------
      if (skipswrad)  then

!---------------------------------------------------------------------
!    calculate shortwave radiative forcing and fluxes using the 
!    exponential-sum-fit parameterization.
!---------------------------------------------------------------------
      else

!---------------------------------------------------------------------

          do indx = 1, size(Sw_output,1)

         !------------------------------------------------------------
         !  determine how aerosols are computed (or not computed)
         !  in the shortwave code
         !------------------------------------------------------------
             if ( .not. do_swaerosol_forcing ) then
              ! when aero forcing is not computed, shortwave is
              ! called once with the value of 'do_swaerosol'
                use_aero = do_swaerosol
             else
              ! when aero forcing is computed, shortwave is called
              ! twice: first time with the value of 'do_swaerosol',
              ! second time with the opposite value
                if (indx .eq. 1) use_aero = do_swaerosol
                if (indx .eq. 2) use_aero = .not.do_swaerosol
             endif

!-----------------------------------------------------------------------
!    
!     sw_cs = .true. turns off the sw cloud radiative effect
!     default; sw_cs = .false.
!-----------------------------------------------------------------------
             camtsw_local= camtsw
        
             if (sw_cs) then
               if (clear_min_sw .EQ. 0 .OR. clear_max_sw .EQ. 0) then
                 camtsw_local= 0.0E+00
               else
                 camtsw_local(:,:,clear_min_sw:clear_max_sw,:) = 0.0E+00
               end if 
             endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!    call swresf with aerosols (if model is being run without) or without
!    aerosols (if model is being run with). save the radiation fluxes 
!    in Sw_output_ad (which does not feed back into the model), but 
!    which may be used to define the aerosol forcing.
!-----------------------------------------------------------------------
             if (use_aero) then
                call swresf (press, pflux, temp, rh2o, deltaz,  &
                             asfc_vis_dir, asfc_nir_dir, &
                             asfc_vis_dif, asfc_nir_dif, &
                             Rad_gases%qo3, Rad_gases%rrvco2,    &   
                             Rad_gases%rrvch4, Rad_gases%rrvn2o, &
                             solflxband, solar_constant_used, &
                             Astro%rrsun, Astro%cosz, Astro%fracday, &
                             camtsw_local, cldsct, cldext, cldasymm, &
                             aeroasymfac, aerosctopdep, aeroextopdep, &
                             Rad_control%do_totcld_forcing, &
                             flag_stoch, Sw_output(indx))
             else
                aerozero = 0.0
                call swresf (press, pflux, temp, rh2o, deltaz,  &
                             asfc_vis_dir, asfc_nir_dir, &
                             asfc_vis_dif, asfc_nir_dif, &
                             Rad_gases%qo3, Rad_gases%rrvco2,    &   
                             Rad_gases%rrvch4, Rad_gases%rrvn2o, &
                             solflxband, solar_constant_used, &
                             Astro%rrsun, Astro%cosz, Astro%fracday, &
                             camtsw_local, cldsct, cldext, cldasymm, &
                             aerozero, aerozero, aerozero, &
                             Rad_control%do_totcld_forcing, &
                             flag_stoch, Sw_output(indx))
             endif

            if (no_sw_cloud_heating_at_levs) then
            !  Sw_output%hsw(:,:,no_heating_min_sw:no_heating_max_sw) = 0.
            Sw_output(indx)%hsw(:,:,no_heating_min_sw:no_heating_max_sw) = Sw_output(indx)%hswcf(:,:,no_heating_min_sw:no_heating_max_sw)
            end if

          enddo  ! indx
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------

      endif  ! skipswrad

!--------------------------------------------------------------------


end subroutine shortwave_driver


!###################################################################
! <SUBROUTINE NAME="shortwave_driver_end">
!  <OVERVIEW>
!   Code that ends shortwave radiation calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine that simply reset shortwave_driver_initialized to false
!  </DESCRIPTION>
!  <TEMPLATE>
!   call shortwave_driver_end
!  </TEMPLATE>
! </SUBROUTINE>
subroutine shortwave_driver_end

!---------------------------------------------------------------------
!    shortwave_driver_end is the destructor for shortwave_driver_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('shortwave_driver_mod',   &
             'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    close out the modules initialized by this module.
!--------------------------------------------------------------------
      call esfsw_driver_end
      call solar_data_driver_end

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!-------------------------------------------------------------------

end subroutine shortwave_driver_end

!#####################################################################

subroutine get_solar_flux_by_band (solarflux, ref)
real,    intent(out)          :: solarflux(:)
logical, intent(in), optional :: ref
!---------------------------------------------------------------------
!    be sure module has been initialized
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('shortwave_driver_mod',   &
             'module has not been initialized', FATAL )
      endif
!---------------------------------------------------------------------
      if (.not. allocated(solflxbandref)) then
        call error_mesg ('shortwave_driver_mod',   &
             'module has not been properly initialized', FATAL )
      endif
!---------------------------------------------------------------------
      if (size(solarflux,1) .ne. size(solflxbandref,1)) then
        call error_mesg ('shortwave_driver_mod',   &
             'incorrect number of bands in get_solar_flux_by_band', FATAL)
      endif
!---------------------------------------------------------------------
      if (present(ref)) then
         if (ref) then
            solarflux = solflxbandref
            return
         endif
      endif

      if (solflxband_initialized) then
         solarflux = solflxband
      else
         call error_mesg ('shortwave_driver_mod',   &
             'time varying solar flux not initialized', FATAL )
      endif
!---------------------------------------------------------------------

end subroutine get_solar_flux_by_band

!#####################################################################

subroutine get_solar_constant (solar_constant)
real,    intent(out) :: solar_constant
!---------------------------------------------------------------------
!    be sure module has been initialized
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('shortwave_driver_mod',   &
             'module has not been initialized', FATAL )
      endif
!---------------------------------------------------------------------
      solar_constant = solar_constant_used
!---------------------------------------------------------------------

end subroutine get_solar_constant
     
!#####################################################################
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine swresf_wrapper ( press, pflux, temp, rh2o, deltaz,  &
                            asfc_vis_dir, asfc_nir_dir, &
                            asfc_vis_dif, asfc_nir_dif, &
                            Rad_gases, Astro, &
                            camtsw, cldsct, cldext, cldasymm, &
                            aeroasymfac, aerosctopdep, aeroextopdep, &
                            do_totcld_forcing, flag_stoch, Sw_output )

real, dimension(:,:,:),     intent(in)    :: press, pflux, temp, rh2o, deltaz
real, dimension(:,:),       intent(in)    :: asfc_vis_dir, asfc_nir_dir, &
                                             asfc_vis_dif, asfc_nir_dif
type(radiative_gases_type), intent(in)    :: Rad_gases   
type(astronomy_type),       intent(in)    :: Astro           
real, dimension(:,:,:,:),   intent(in)    :: camtsw
real, dimension(:,:,:,:,:), intent(in)    :: cldsct, cldext, cldasymm
real,dimension(:,:,:,:),    intent(in)    :: aeroasymfac, aerosctopdep, aeroextopdep
logical,                    intent(in)    :: do_totcld_forcing
integer,                    intent(in)    :: flag_stoch
type(sw_output_type),       intent(inout) :: Sw_output

real, dimension(size(Astro%cosz,1), &
                size(Astro%cosz,2)) :: cosz, fracday
!real    :: camtsw_local(size(camtsw,1),size(camtsw,2),size(camtsw,3),size(camtsw,4))
integer :: kx

!%    camtsw_local= camtsw
!%
!%! sw_cs = .true. turns off the sw cloud radiative effect
!%! default; sw_cs = .false.
!%    if (sw_cs) then
!%      if (clear_min_sw .EQ. 0 .OR. clear_max_sw .EQ. 0) then
!%        camtsw_local= 0.0E+00
!%      else
!%        camtsw_local(:,:,clear_min_sw:clear_max_sw,:) = 0.0E+00
!%      end if 
!%    endif

   kx = size (press,3)

   call swresf (press (:,:,1:kx),   &
                pflux (:,:,1:kx+1), &
                temp  (:,:,1:kx),   &
                rh2o  (:,:,1:kx),   &
                deltaz(:,:,1:kx),   &
                asfc_vis_dir, asfc_nir_dir, &
                asfc_vis_dif, asfc_nir_dif, &
                Rad_gases%qo3, Rad_gases%rrvco2,    &
                Rad_gases%rrvch4, Rad_gases%rrvn2o, &
                solflxband, solar_constant_used, &
                Astro%rrsun, Astro%cosz, Astro%fracday, &
                camtsw, cldsct, cldext, cldasymm, &
                aeroasymfac, aerosctopdep, aeroextopdep, &
                do_totcld_forcing, flag_stoch, Sw_output)

!-------------------------------------------------------------------------------
 ! mtp: replace sw-all-sky heating rates by clear sky heating rates from
 ! level no_cloud_min to level no_cloud_max
 ! level 1=top layer, level 32 = surface layer (for 32 layer atmosphere)
 ! levl 20 = 663 hPa

!if (no_sw_cloud_heating_at_levs) then
!!  Sw_output%hsw(:,:,no_heating_min_sw:no_heating_max_sw) = 0.
!Sw_output%hsw(:,:,no_heating_min_sw:no_heating_max_sw) = Sw_output%hswcf(:,:,no_heating_min_sw:no_heating_max_sw)
!end if

end subroutine swresf_wrapper

!####################################################################

subroutine compute_time_offset ( Time, entry_date, &
                                 Entry_time, Offset, negative_offset )

type(time_type),       intent(in)  :: Time
integer, dimension(6), intent(in)  :: entry_date
type(time_type),       intent(out) :: Entry_time, Offset
logical,               intent(out) :: negative_offset

!--------------------------------------------------------------------
!  compute the time offset (as a time_type variable) between
!  date and Time (e.g. date - Time). The difference will always
!  be positive. If date < Time then negative_offset = true.
!--------------------------------------------------------------------

      if (entry_date(1) == 1 .and. &
          entry_date(2) == 1 .and. &
          entry_date(3) == 1 .and. &
          entry_date(4) == 0 .and. &
          entry_date(5) == 0 .and. &
          entry_date(6) == 0 ) then
             Entry_time = Time
      else
             Entry_time = set_date (entry_date(1), &
                                    entry_date(2), &
                                    entry_date(3), &
                                    entry_date(4), &
                                    entry_date(5), &
                                    entry_date(6))
      endif

      Offset = Entry_time - Time
      if (Time > Entry_time) then
         negative_offset = .true.
      else
         negative_offset = .false.
      endif

end subroutine compute_time_offset

!####################################################################


                end module shortwave_driver_mod

