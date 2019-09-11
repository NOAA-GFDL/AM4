                       module rad_output_file_mod

! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Module that provides subroutines to write radiation output to
!  history file
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

!   shared modules:

use mpp_mod,           only: input_nml_file
use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, NOTE, close_file, string, lowercase
use time_manager_mod,  only: time_manager_init, time_type, operator(>)
use diag_manager_mod,  only: register_diag_field, diag_manager_init, &
                             send_data, get_diag_field_id, &
                             DIAG_FIELD_NOT_FOUND
use diag_axis_mod,     only: get_axis_num
use diag_data_mod,     only: CMOR_MISSING_VALUE
use constants_mod,     only: constants_init, GRAV, WTMAIR, WTMOZONE, pi

use aerosol_types_mod, only:  aerosol_type

use atmos_cmip_diag_mod, only: register_cmip_diag_field_3d, &
                               register_cmip_diag_field_2d, &
                               send_cmip_data_3d, &
                               cmip_diag_id_type, &
                               query_cmip_diag_id

!  radiation package shared modules:

use radiation_driver_types_mod, only:  radiation_control_type, &
                                       rad_output_type

use aerosolrad_types_mod, only:  aerosolrad_control_type, &
                                 aerosolrad_diag_type

use cloudrad_types_mod,   only:  cld_specification_type

use esfsw_driver_mod,     only:  esfsw_bands
use shortwave_types_mod,  only:  sw_output_type
use longwave_types_mod,   only:  lw_output_type

!--------------------------------------------------------------------

implicit none
private

!------------------------------------------------------------------
!    rad_output_file_mod writes an output file containing an assort-
!    ment of variables related to the sea_esf_rad radiation package.
!    this is an optionally-generated file, which may be used to sup-
!    plement the standard diagnostic model output files. NOTE THAT
!    THIS FILE IS GENERATED ONLY ON RADIATION TIMESTEPS, SO THAT WHEN 
!    SW FLUXES ARE BEING RENORMALIZED ON EACH PHYSICS STEP, VARIABLES 
!    IN THIS FILE RELATED TO SW RADIATION WILL NOT REFLECT THE EFFECTS 
!    OF THE RENORMALIZATION.
!------------------------------------------------------------------


!-------------------------------------------------------------------
!----------- version number for this module ------------------------

character(len=128)  :: version = &
'$Id$'
character(len=128)  :: tagname =  '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public   &
         rad_output_file_init, write_rad_output_file,   &
         rad_output_file_end

private   &

!   called from rad_output_file_init
        register_fields


!-------------------------------------------------------------------
!-------- namelist  ---------

logical :: write_data_file=.false.  ! data file to be written  ?


namelist  / rad_output_file_nml /  &
                                  write_data_file
     
!------------------------------------------------------------------
!----public data------


!------------------------------------------------------------------
!----private data------

!--------------------------------------------------------------------
!    DU_factor is Dobson units per (kg/m2). the first term is 
!    (avogadro's number/loeschmidt's number at STP = vol./kmol of an 
!    ideal gas at STP). second term = mol wt o3 . third term is units 
!    conversion. values are chosen from US Standard Atmospheres, 1976.
!--------------------------------------------------------------------
real, parameter  :: DU_factor =    &
                            (6.022169e26/2.68684e25)/(47.9982)*1.0e5
real, parameter  :: DU_factor2 = DU_factor/GRAV
                                   ! Dobson units per (kg/kg * dyn/cm2) 
                                   ! Dobson units per (kg/kg * N  /m2) 

!--------------------------------------------------------------------
! netcdf diagnostics field variables
!---------------------------------------------------------------------
character(len=16), parameter       :: mod_name='radiation'
real                               :: missing_value = -999.
integer, dimension(:), allocatable :: id_aerosol, id_aerosol_column
!integer, dimension(:), allocatable :: id_aerosol, id_aerosol_column, &
!                                     id_absopdep, id_absopdep_column, &
integer, dimension(:,:), allocatable :: id_absopdep,  &
                                        id_absopdep_column, &
                                      id_extopdep, id_extopdep_column
integer, dimension(:,:), allocatable :: id_asymdep, id_asymdep_column
integer, dimension(2)              :: id_lw_absopdep_vlcno_column, &
                                      id_lw_extopdep_vlcno_column, &
                                      id_lwext_vlcno, id_lwssa_vlcno, &
                                      id_lwasy_vlcno, id_lw_xcoeff_vlcno
integer, dimension(2)              :: id_absopdep_vlcno_column, &
                                      id_extopdep_vlcno_column, &
                                      id_swext_vlcno, id_swssa_vlcno, &
                                      id_swasy_vlcno, id_sw_xcoeff_vlcno
! change by Xianglei Huang
integer, dimension(7)              :: id_lw_bdyflx_clr, id_lw_bdyflx
integer                            :: id_bT_Aqua31
! end of change
integer, dimension(4)              :: id_sw_bdyflx_clr, id_sw_bdyflx
integer                            :: id_swheat_vlcno
integer                            :: id_sulfate_col_cmip,  &
                                      id_sulfate_cmip
integer, dimension(:), allocatable :: id_aerosol_fam, &
                                      id_aerosol_fam_column    
integer, dimension(:,:), allocatable :: id_absopdep_fam,  &
                                        id_absopdep_fam_column, &
                                        id_extopdep_fam,  &
                                        id_extopdep_fam_column
integer, dimension(:,:), allocatable :: id_asymdep_fam,  &
                                        id_asymdep_fam_column
integer                            :: id_radswp, id_radp, id_temp, &
                                      id_rh2o, id_qo3, id_qo3_col,  &
                                      id_qo3v, &
                                      id_cmxolw, id_crndlw, id_flxnet, &
                                      id_fsw, id_ufsw, id_psj, &
                                      id_dfsw, &
                                      id_tmpsfc, id_cvisrfgd_dir,  &
                                      id_cirrfgd_dir, &
                                      id_cvisrfgd_dif, id_cirrfgd_dif, &
                                      id_radswpcf,  &
                                      id_cldwater, id_cldice,  &
                                      id_cldarea, &
                                      id_radpcf, id_flxnetcf, &
                                      id_fswcf, id_ufswcf, id_pressm,  &
                                      id_dfswcf, &
                                      id_phalfm, id_pfluxm, &
                                      id_dphalf, id_dpflux, &
                                      id_ptop

type(cmip_diag_id_type)  :: ID_o3, ID_ec550aer, ID_concso4, ID_concsoa
integer                  :: id_loadso4, id_sconcso4, id_loadsoa, id_sconcsoa, &
                            id_od550aer, id_od550lt1aer, id_abs550aer, id_od870aer

! cmip names for select tracer families
! also partial long_names and standard_names
integer, dimension(6) :: cmip_family_mapping
integer, dimension(6) :: id_cmipload, id_cmipsconc
type(cmip_diag_id_type), dimension(6) :: ID_cmipconc
character(len=8), dimension(6) :: cmip_names = (/"oa  ","poa ","soa ","bc  ","dust","ss  "/)
character(len=64), dimension(6) :: cmip_longnames = &
                                  (/"Dry Aerosol Organic Matter          ", &
                                    "Dry Aerosol Primary Organic Matter  ", &
                                    "Dry Aerosol Secondary Organic Matter", &
                                    "Black Carbon Aerosol                ", &
                                    "Dust                                ", &
                                    "Seasalt                             "/)
character(len=64), dimension(6) :: cmip_stdnames = &
                                  (/"particulate_organic_matter          ", &
                                    "primary_particulate_organic_matter  ", &
                                    "secondary_particulate_organic_matter", &
                                    "black_carbon                        ", &
                                    "dust                                ", &
                                    "seasalt                             "/)

!---------------------------------------------------------------------
!    miscellaneous variables
!---------------------------------------------------------------------
integer :: nso4, nsoa, naero, npm25, nvis, n870
integer :: naerosol=0                      ! number of active aerosols
logical :: module_is_initialized= .false.  ! module initialized ?
integer, parameter              :: N_DIAG_BANDS = 10
character(len=16), dimension(N_DIAG_BANDS) ::   &
                     band_suffix = (/ '_vis', '_nir', '_con',  &
                                      '_bd5', '_bd6', '_870', &
! +++ pag 11/13/2009
                                      '_340','_380','_440','_670' /)
! Properties at specific wavelength are in fact averaged over a band which
! are defined in the input file:
! /home/pag/fms/radiation/esf_sw_input_data_n38b18_1992_version_ckd2.1.lean.nov89.ref
!
! 309 < 340 < 364
! 364 < 380 < 406
! 406 < 440 < 448
! 500 < vis < 600
! 600 < 670 < 685
! 685 < 870 < 870
! 870 < nir < 1219
! In the file, each band is defined by its end wavelength in wavenumber.
! For 340nm, the endband is 1.e4/0.309=32400 (11th wavenumber among the 18
! specified in the input file)
! +++ pag 11/13/2009

! shortwave band indices
!   w550_band_indx = visible
!   one_micron_indx = near infra red

integer :: w550_band_indx, one_micron_indx


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!#####################################################################
! <SUBROUTINE NAME="rad_output_file_init">
!  <OVERVIEW>
!   Constructor of rad_output_file module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to initialize and set up rad_output_file module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  rad_output_file_init (axes, Time, names)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="integer">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
!  <IN NAME="names" TYPE="character">
!   aerosol names
!  </IN>
! </SUBROUTINE>
!
subroutine rad_output_file_init (axes, Time, names, family_names, &
                                 Rad_control, Aerosolrad_control)

!--------------------------------------------------------------------
!    rad_output_file_init is the constructor for rad_output_file_mod.
!--------------------------------------------------------------------

integer, dimension(4),          intent(in)    :: axes
type(time_type),                intent(in)    :: Time
character(len=*), dimension(:), intent(in)    :: names
character(len=*), dimension(:), intent(in)    :: family_names
type(radiation_control_type),   intent(in)    :: Rad_control
type(aerosolrad_control_type),  intent(in)    :: Aerosolrad_control

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    these variables are present when running the gcm, not present
!    when running the standalone code.
!  
!       axes      diagnostic variable axes for netcdf files
!       Time      current time [ time_type(days, seconds) ]
!       names     names of active aerosols
!       family_names 
!                 names of active aerosol families
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer   :: unit, io, ierr, logunit
      integer   :: nfields

!---------------------------------------------------------------------
!   local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!        nfields         number of active aerosol fields
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
      call constants_init
      call diag_manager_init  
      call time_manager_init 

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=rad_output_file_nml, iostat=io)
      ierr = check_nml_error(io,'rad_output_file_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=rad_output_file_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'rad_output_file_nml')
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
                         write (logunit, nml=rad_output_file_nml)

!--------------------------------------------------------------------
!    if running gcm, continue on if data file is to be written. 
!--------------------------------------------------------------------
        if (write_data_file) then

!---------------------------------------------------------------------
!    register the diagnostic fields for output.
!---------------------------------------------------------------------
          nfields = size(names(:))
          call register_fields (Time, axes, nfields, names, family_names, &
                                Rad_control%do_totcld_forcing, &
                                Aerosolrad_control%volcanic_lw_aerosols, &
                                Aerosolrad_control%volcanic_sw_aerosols)

!--------------------------------------------------------------------
!     determine index for visible and infrared short bands
!--------------------------------------------------------------------
          call esfsw_bands ( w550_band_indx=w550_band_indx, &
                             one_micron_indx=one_micron_indx )

        endif

!--------------------------------------------------------------------
!    mark the module as initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------


end subroutine rad_output_file_init



!################################################################
! <SUBROUTINE NAME="write_rad_output_file">
!  <OVERVIEW>
!   write_rad_output_file produces a netcdf output file containing
!   the user-specified radiation-related variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!   write_rad_output_file produces a netcdf output file containing
!   the user-specified radiation-related variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call write_rad_output_file (is, ie, js, je, Atmos_input, Surface, &
!                                  Rad_output, &
!                                  Sw_output, Lw_output, Rad_gases, &
!                                  Cld_spec,  &
!                                  Time_diag, Time, aerosol_in)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending subdomain i,j indices of data 
!   in the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable containing atmos-
!                        pheric input data for the radiation package 
!                        on the model grid
!  </IN>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface input data to radiation package
!  </IN>
!  <IN NAME="Rad_output" TYPE="rad_output_type">
!   rad_output_type variable containing radiation
!                        output data needed by other modules
!  </IN>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   sw_output_type variable containing shortwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid
!  </IN>
!  <IN NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid 
!  </IN>
!  <IN NAME="Cldrad_pros" TYPE="cldrad_properties_type">
!   cldrad_properties_type variable containing 
!                        cloud radiative property input data for the 
!                        radiation package on the model grid
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_diagnostics_type">
!   cld_specification_type variable containing cloud
!                        microphysical data
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!   time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ]
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
!  <IN NAME="aerosol_in" TYPE="real">
!   optional aerosol data
!  </IN>
! </SUBROUTINE>
!
subroutine write_rad_output_file (is, ie, js, je, tmpsfc, &
                                  asfc_nir_dir, asfc_vis_dir, &
                                  asfc_nir_dif, asfc_vis_dif, &
                                  press, pflux, phalf, &
                                  temp, rh2o, qo3, deltaz,  &
                                  Rad_control, Aerosolrad_control, &
                                  Rad_output, Sw_output, Lw_output,    &
                                  Cld_spec, Time_diag, Time, &
                                  Aerosol, Aerosolrad_diags)

!----------------------------------------------------------------
!    write_rad_output_file produces a netcdf output file containing
!    the user-specified radiation-related variables.
!----------------------------------------------------------------

integer,                      intent(in)            ::  is, ie, js, je
real, dimension(:,:),         intent(in)            ::  tmpsfc, &
                                                        asfc_nir_dir, &
                                                        asfc_vis_dir, &
                                                        asfc_nir_dif, &
                                                        asfc_vis_dif
real, dimension(:,:,:),       intent(in)            ::  press, pflux, &
                                                        phalf, temp, &
                                                        rh2o, qo3, deltaz
type(radiation_control_type),  intent(in)           ::  Rad_control
type(aerosolrad_control_type), intent(in)           ::  Aerosolrad_control
type(rad_output_type),        intent(in)            ::  Rad_output
type(sw_output_type),         intent(in)            ::  Sw_output
type(lw_output_type),         intent(in)            ::  Lw_output
type(cld_specification_type), intent(in)            ::  Cld_spec      
type(time_type),              intent(in)            ::  Time_diag, Time
type(aerosol_type),           intent(in), optional  ::  Aerosol
type(aerosolrad_diag_type),   intent(in), optional  ::  Aerosolrad_diags

!------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      Atmos_input       atmos_input_type variable containing atmos-
!                        pheric input data for the radiation package 
!                        on the model grid
!      Surface           surface input fields to radiation package
!                        [ surface_type ]
!      qo3               ozone mixing ratio [ g / g ]
!      Rad_output        rad_output_type variable containing radiation
!                        output data needed by other modules
!      Sw_output         sw_output_type variable containing shortwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid
!      Lw_output         lw_output_type variable containing longwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid 
!      Cld_spec          cld_specification_type variable containing 
!                        cloud specification input data for the 
!                        radiation package on the model grid
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ]  
!      Time              current time [ time_type(days, seconds) ]

!
!  intent(in), optional variables:
!
!      aerosol_in        active aerosol distributions
!                        [ kg / m**2 ]
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      real, dimension (size(press,1),   &
                       size(press,2) )  ::  &
                           psj, qo3_col, ptop

      real, dimension (size(phalf,1),    &
                       size(phalf,2), &
                       size(phalf,3)  ) ::   &
                          fsw, ufsw, fswcf, ufswcf, flxnet, flxnetcf

      real, dimension (size(press,1), &
                       size(press,2), &
                       size(press,3) ) ::  &
                       cmxolw, crndlw, radp, radswp, v_heat, &
                       radpcf, radswpcf, dphalf, dpflux

      real, dimension (size(press,1),    &
                       size(press,2), 4)  :: &
                                bdy_flx_mean, bdy_flx_clr_mean

      real, dimension(size(press,1),    &
                       size(press,2)) :: bT

      real, dimension(:,:,:),   allocatable :: aerosol_col
      real, dimension(:,:,:,:),   allocatable :: extopdep_col
      real, dimension(:,:,:,:),   allocatable :: absopdep_col
      real, dimension(:,:,:),   allocatable :: lw_extopdep_vlcno_col
      real, dimension(:,:,:),   allocatable :: lw_absopdep_vlcno_col
      real, dimension(:,:,:),   allocatable :: extopdep_vlcno_col
      real, dimension(:,:,:),   allocatable :: absopdep_vlcno_col
      real, dimension(:,:,:,:),   allocatable :: absopdep_fam_col
      real, dimension(:,:,:,:),   allocatable :: extopdep_fam_col
      real, dimension(:,:,:),   allocatable :: aerosol_fam_col
      real, dimension(:,:,:,:,:), allocatable :: absopdep_fam
      real, dimension(:,:,:,:,:), allocatable :: extopdep_fam
      real, dimension(:,:,:,:), allocatable :: aerosol_fam
      real, dimension(:,:,:,:),   allocatable :: asymdep_col
      real, dimension(:,:,:,:,:), allocatable :: asymdep_fam
      real, dimension(:,:,:,:),   allocatable :: asymdep_fam_col
      real, dimension(:,:,:,:), allocatable :: sum1
      real, dimension(:,:,:), allocatable :: sum2

      logical   :: used, Lasymdep
      integer   :: kerad ! number of model layers
      integer   :: n, k, na, nfamilies, nl
      integer   :: nv
      integer   :: co_indx, bnd_indx
      integer   :: ncmip

!----------------------------------------------------------------------
!  local variables:
!
!      tmpsfc         surface temperature [ deg K ]
!      psj            surface pressure [ hPa ]
!      cvisrfgd       surface visible light albedo [ dimensionless ]
!      cirrfgd        surface ir albedo [ dimensionless ]
!      tot_clds       total column isccp clouds [ percent ]
!      cld_isccp_hi   number of isccp high clouds [ percent ]
!      cld_isccp_mid  number of isccp middle clouds [ percent ]
!      cld_isccp_low  number of isccp low clouds [ percent ]
!      qo3_col        ozone column [ DU ]
!      fsw            net shortwave flux [ W / m**2 ]
!      ufsw           upward shortwave flux [ W / m**2 ]
!      fswcf          net sw flux in the absence of clouds [ W / m**2 ]
!      ufswcf         upward sw flux in absence of clouds [ W / m**2]
!      flxnet         net longwave flux [ W / m**2 ]
!      flxnetcf       net lw flux in the absence of clouds [ W / m**2 ]
!      phalfm         model interface level pressure [ Pa ]
!      pfluxm         avg of adjacent model level pressures [ Pa ]
!      temp           temperature [ deg K ]
!      rh2o           water vapor specific humidity [ g / g ]
!      heatra         lw heating rate [ deg K / day ]
!      heatracf       lw heating rate without cloud [ deg K / day ]
!      cmxolw         amount of maximal overlap clouds [ percent]
!      crndlw         amount of ramndom overlap clouds [ percent]
!      radp           lw + sw heating rate [ deg K / sec ]
!      radswp         sw heating rate [ deg K / sec ]
!      radpcf         lw + sw heating rate w/o clouds [ deg K / sec ]
!      radswpcf       sw heating rate w/o clouds [ deg K / sec ]
!      pressm         pressure at model levels [ Pa ]
!      aerosol_col
!      used
!      kerad
!      n,k
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_output_file_mod', &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    if the file is not to be written, do nothing.
!--------------------------------------------------------------------
      if (write_data_file) then
      if (Time_diag > Time) then

        Lasymdep = .false.
        if (any(id_asymdep_fam(:,:)  > 0) .or. &
            any(id_asymdep_fam_column(:,:)  > 0 )) Lasymdep = .true.
            
!--------------------------------------------------------------------
!    if the file is to be written, define the number of model layers.
!--------------------------------------------------------------------
        kerad = ubound(temp,3)

!--------------------------------------------------------------------
!    retrieve the desired fields from the input derived data types.
!--------------------------------------------------------------------
        psj   (:,:)     = phalf(:,:,kerad+1)

        ! for reproducibility
        radp   = Rad_output%tdt_rad(is:ie,js:je,:)/float(1)
        radswp = Rad_output%tdtsw  (is:ie,js:je,:)/float(1)
        fsw  = Sw_output%fsw (:,:,:)/float(1)
        ufsw = Sw_output%ufsw(:,:,:)/float(1)
        bdy_flx_mean = Sw_output%bdy_flx(:,:,:)/float(1)

!BW     cirrfgd_dir(:,:)    = Surface%asfc_nir_dir(:,:)
!BW     cvisrfgd_dir(:,:)   = Surface%asfc_vis_dir(:,:)
!BW     cirrfgd_dif(:,:)    = Surface%asfc_nir_dif(:,:)
!BW     cvisrfgd_dif(:,:)   = Surface%asfc_vis_dif(:,:)
        flxnet(:,:,:)   = Lw_output%flxnet(:,:,:)
!BW     qo3(:,:,:)      = Rad_gases%qo3(:,:,:)
        cmxolw(:,:,:)   = 100.0*Cld_spec%cmxolw(:,:,:)
        crndlw(:,:,:)   = 100.0*Cld_spec%crndlw(:,:,:)

          do k = 1,kerad
            dphalf(:,:,k)   = phalf(:,:,k+1) - phalf(:,:,k)
            dpflux(:,:,k)   = pflux(:,:,k+1) - pflux(:,:,k)
          enddo
          ptop(:,:) = 0.01*phalf(:,:,1)

        if (Rad_control%do_totcld_forcing) then 
          ! for reproducibility
          fswcf  = Sw_output%fswcf(:,:,:)/float(1)
          ufswcf = Sw_output%ufswcf(:,:,:)/float(1)
          radpcf = Rad_output%tdt_rad_clr(is:ie,js:je,:)/float(1)
          radswpcf = Rad_output%tdtsw_clr  (is:ie,js:je,:)/float(1)
          bdy_flx_clr_mean = Sw_output%bdy_flx_clr(:,:,:)/float(1)
          flxnetcf(:,:,:) = Lw_output%flxnetcf(:,:,:)
        endif

!---------------------------------------------------------------------
!    calculate the column ozone in DU (Dobson units). convert from 
!    (kg/kg) * (N/m2) to DU (1DU = 2.687E16 molec cm^-2).
!---------------------------------------------------------------------
        qo3_col(:,:) = 0.
        do k = 1,size(qo3,3)
          qo3_col(:,:) = qo3_col(:,:) + &
                         qo3(:,:,k)*(pflux(:,:,k+1) - pflux(:,:,k))
        end do
        qo3_col(:,:) = qo3_col(:,:)*DU_factor2

!---------------------------------------------------------------------
!    define the aerosol fields and calculate the column aerosol. 
!---------------------------------------------------------------------
        if (Aerosolrad_control%do_aerosol) then
          allocate ( aerosol_col(size(Aerosol%aerosol, 1), &
                                 size(Aerosol%aerosol, 2), &
                                 size(Aerosol%aerosol, 4)) )       
          aerosol_col(:,:,:) = SUM (Aerosol%aerosol(:,:,:,:), 3)
!         if (Aerosolrad_control%do_swaerosol) then
            allocate ( extopdep_col(size(Aerosolrad_diags%extopdep  , 1), &
                                    size(Aerosolrad_diags%extopdep  , 2), &
                               size(Aerosolrad_diags%extopdep  , 4), &
                                    N_DIAG_BANDS) )
            extopdep_col(:,:,:,:) =    &
                          SUM (Aerosolrad_diags%extopdep  (:,:,:,:,:), 3)
            allocate ( absopdep_col(size(Aerosolrad_diags%absopdep  , 1), &
                                    size(Aerosolrad_diags%absopdep  , 2), &
                              size(Aerosolrad_diags%absopdep  , 4), &
                                   N_DIAG_BANDS) )
            absopdep_col(:,:,:,:) =    &
                           SUM (Aerosolrad_diags%absopdep  (:,:,:,:,:), 3)
            if ( Lasymdep ) then 
              allocate ( asymdep_col(size(Aerosolrad_diags%asymdep  , 1), &
                                 size(Aerosolrad_diags%asymdep  ,     2), &
                                 size(Aerosolrad_diags%asymdep  ,     4), &
                                      N_DIAG_BANDS) )
              asymdep_col(:,:,:,:) =  SUM ( &
                            Aerosolrad_diags%asymdep(:,:,:,:,:) &
                                *(Aerosolrad_diags%extopdep(:,:,:,:,:)-  &
                                    Aerosolrad_diags%absopdep(:,:,:,:,:)), 3) &
                             /(1.e-30+SUM(Aerosolrad_diags%extopdep(:,:,:,:,:)&
                                     -Aerosolrad_diags%absopdep(:,:,:,:,:), 3))
            endif
            if (Aerosolrad_control%volcanic_sw_aerosols) then
              allocate ( extopdep_vlcno_col(   &
                           size(Aerosolrad_diags%extopdep_vlcno  , 1), &
                           size(Aerosolrad_diags%extopdep_vlcno  , 2),3))
              extopdep_vlcno_col(:,:,:) =    &
                        SUM (Aerosolrad_diags%extopdep_vlcno  (:,:,:,:), 3)
              allocate ( absopdep_vlcno_col(    &
                            size(Aerosolrad_diags%absopdep_vlcno  , 1), &
                            size(Aerosolrad_diags%absopdep_vlcno  , 2),3))
              absopdep_vlcno_col(:,:,:) =    &
                        SUM (Aerosolrad_diags%absopdep_vlcno  (:,:,:,:), 3)
            endif
              
            if (Aerosolrad_control%volcanic_lw_aerosols) then
              allocate ( lw_extopdep_vlcno_col(   &
                         size(Aerosolrad_diags%lw_extopdep_vlcno  , 1), &
                         size(Aerosolrad_diags%lw_extopdep_vlcno  , 2),2))
              lw_extopdep_vlcno_col(:,:,:) =    &
                    SUM (Aerosolrad_diags%lw_extopdep_vlcno  (:,:,:,:), 3)
              allocate ( lw_absopdep_vlcno_col(    &
                          size(Aerosolrad_diags%lw_absopdep_vlcno  , 1), &
                          size(Aerosolrad_diags%lw_absopdep_vlcno  , 2),2))
              lw_absopdep_vlcno_col(:,:,:) =    &
                    SUM (Aerosolrad_diags%lw_absopdep_vlcno  (:,:,:,:), 3)
            endif
        endif
        
!---------------------------------------------------------------------
!    define the aerosol family output fields.
!---------------------------------------------------------------------
        if (Aerosolrad_control%do_aerosol) then
          nfamilies = size(Aerosol%family_members,2)
          if (nfamilies > 0) then
            allocate (aerosol_fam (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    size(Aerosol%aerosol,3), &
                                    nfamilies))
            allocate (aerosol_fam_col (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    nfamilies))
            allocate (extopdep_fam (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    size(Aerosol%aerosol,3), &
                                    nfamilies, N_DIAG_BANDS))
            allocate (absopdep_fam (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    size(Aerosol%aerosol,3), &
                                    nfamilies, N_DIAG_BANDS))
            allocate (extopdep_fam_col (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    nfamilies, N_DIAG_BANDS))
            allocate (absopdep_fam_col (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    nfamilies, N_DIAG_BANDS))
            aerosol_fam = 0.
            aerosol_fam_col = 0.
            extopdep_fam = 0.
            absopdep_fam = 0.
            extopdep_fam_col = 0.
            absopdep_fam_col = 0.

            if ( Lasymdep ) then
              allocate (asymdep_fam ( size(Aerosol%aerosol,1), &
                                      size(Aerosol%aerosol,2), &
                                      size(Aerosol%aerosol,3), &
                                      nfamilies, N_DIAG_BANDS))
              allocate (asymdep_fam_col ( size(Aerosol%aerosol,1), &
                                          size(Aerosol%aerosol,2), &
                                          nfamilies, N_DIAG_BANDS))
              allocate (sum1 (       size(Aerosol%aerosol,1), &
                                     size(Aerosol%aerosol,2), &
                                     size(Aerosol%aerosol,3), &
                                     N_DIAG_BANDS))
              allocate (sum2 (       size(Aerosol%aerosol,1), &
                                     size(Aerosol%aerosol,2), &
                                     N_DIAG_BANDS))
              asymdep_fam = 0.
              asymdep_fam_col = 0.
            endif

            do n = 1, nfamilies                      
              do na = 1, naerosol                
                if (Aerosol%family_members(na,n)) then
                  aerosol_fam(:,:,:,n) = aerosol_fam(:,:,:,n) +  &
                                         Aerosol%aerosol(:,:,:,na)
                  aerosol_fam_col(:,:,n) = aerosol_fam_col(:,:,n) +  &
                                         aerosol_col(:,:,na)
                  do nl = 1,N_DIAG_BANDS
                    extopdep_fam(:,:,:,n,nl) = extopdep_fam(:,:,:,n,nl) +  &
                                      Aerosolrad_diags%extopdep(:,:,:,na,nl)
                    extopdep_fam_col(:,:,n,nl) = extopdep_fam_col(:,:,n,nl) +  &
                                      extopdep_col(:,:,na,nl)
                    absopdep_fam(:,:,:,n,nl) = absopdep_fam(:,:,:,n,nl) +  &
                                      Aerosolrad_diags%absopdep(:,:,:,na,nl)
                    absopdep_fam_col(:,:,n,nl) = absopdep_fam_col(:,:,n,nl) +  &
                                      absopdep_col(:,:,na,nl)
                  end do ! (nl)
                endif
              end do ! (na)
              if ( Lasymdep ) then
                sum1(:,:,:,:)=1.e-30
                sum2(:,:,:)=1.e-30
                do na = 1, naerosol                
                  if (Aerosol%family_members(na,n)) then
                    do nl = 1,N_DIAG_BANDS
                      asymdep_fam(:,:,:,n,nl) = asymdep_fam(:,:,:,n,nl) &
                          + Aerosolrad_diags%asymdep(:,:,:,na,nl) &
                          *( Aerosolrad_diags%extopdep(:,:,:,na,nl) &
                          -Aerosolrad_diags%absopdep(:,:,:,na,nl))
                      sum1(:,:,:,nl) = sum1(:,:,:,nl) &
                          + Aerosolrad_diags%extopdep(:,:,:,na,nl) &
                          - Aerosolrad_diags%absopdep(:,:,:,na,nl)
                      asymdep_fam_col(:,:,n,nl) = asymdep_fam_col(:,:,n,nl) &
                          + asymdep_col(:,:,na,nl)&
                          *(extopdep_col(:,:,na,nl)-absopdep_col(:,:,na,nl))
                      sum2(:,:,nl) = sum2(:,:,nl)&
                          + extopdep_col(:,:,na,nl) - absopdep_col(:,:,na,nl)
                    end do ! (nl)
                  endif
                end do ! (na)

                asymdep_fam = max (0., min(1., asymdep_fam))
                sum1 = max(1.e-30,min(1.,sum1))
                asymdep_fam_col = max (0., min(1., asymdep_fam_col))
                sum2 = max (1.e-30, min(1.,sum2))
                do nl = 1,N_DIAG_BANDS
                  asymdep_fam(:,:,:,n,nl) = asymdep_fam(:,:,:,n,nl)/sum1(:,:,:,nl)
                  asymdep_fam_col(:,:,n,nl)=asymdep_fam_col(:,:,n,nl)/sum2(:,:,nl)
                enddo
              endif

              if (Aerosol%family_members(naerosol+1,n)) then
                if (Aerosolrad_control%volcanic_sw_aerosols) then
                  extopdep_fam_col(:,:,n,1) = extopdep_fam_col(:,:,n,1) +  &
                                              extopdep_vlcno_col(:,:,1)
                  absopdep_fam_col(:,:,n,1) = absopdep_fam_col(:,:,n,1) +  &
                                              absopdep_vlcno_col(:,:,1)
                  extopdep_fam_col(:,:,n,2) = extopdep_fam_col(:,:,n,2) +  &
                                              extopdep_vlcno_col(:,:,2)
                  absopdep_fam_col(:,:,n,2) = absopdep_fam_col(:,:,n,2) +  &
                                              absopdep_vlcno_col(:,:,2)
                  extopdep_fam_col(:,:,n,6) = extopdep_fam_col(:,:,n,6) +  &
                                              extopdep_vlcno_col(:,:,3)
                  absopdep_fam_col(:,:,n,6) = absopdep_fam_col(:,:,n,6) +  &
                                              absopdep_vlcno_col(:,:,3)
                endif
                if (Aerosolrad_control%volcanic_lw_aerosols) then
                  extopdep_fam_col(:,:,n,4) = extopdep_fam_col(:,:,n,4) +  &
                                              lw_extopdep_vlcno_col(:,:,1)
                  absopdep_fam_col(:,:,n,4) = absopdep_fam_col(:,:,n,4) +  &
                                              lw_absopdep_vlcno_col(:,:,1)
                  extopdep_fam_col(:,:,n,5) = extopdep_fam_col(:,:,n,5) +  &
                                              lw_extopdep_vlcno_col(:,:,2)
                  absopdep_fam_col(:,:,n,5) = absopdep_fam_col(:,:,n,5) +  &
                                              lw_absopdep_vlcno_col(:,:,2)
                endif
              endif
            enddo ! (n)

            do n = 1,nfamilies
              if (id_aerosol_fam(n)  > 0 ) then
                used = send_data (id_aerosol_fam(n),  &
                                  aerosol_fam(:,:,:,n)/deltaz(:,:,:),   &
                                  Time_diag, is, js, 1)
              endif
              if (id_aerosol_fam_column(n)  > 0 ) then
                used = send_data (id_aerosol_fam_column(n),     &
                                  aerosol_fam_col(:,:,n), Time_diag, is, js)
              endif

              do nl=1,N_DIAG_BANDS
                if (id_extopdep_fam(n,nl)  > 0 ) then
                  used = send_data (id_extopdep_fam(n,nl),    &
                                    extopdep_fam  (:,:,:,n,nl), &
                                    Time_diag, is, js, 1)
                endif
                if (id_extopdep_fam_column(n,nl)  > 0 ) then
                  used = send_data (id_extopdep_fam_column(n,nl),     &
                                   extopdep_fam_col(:,:,n,nl), Time_diag, is, js)
                endif
                if (id_absopdep_fam(n,nl)  > 0 ) then
                  used = send_data (id_absopdep_fam(n,nl),    &
                                    absopdep_fam  (:,:,:,n,nl), &
                                    Time_diag, is, js, 1)
                endif
                if (id_absopdep_fam_column(n,nl)  > 0 ) then
                  used = send_data (id_absopdep_fam_column(n,nl),     &
                               absopdep_fam_col(:,:,n,nl), Time_diag, is, js)
                endif
                if (id_asymdep_fam(n,nl)  > 0 ) then
                  used = send_data (id_asymdep_fam(n,nl),    &
                                    asymdep_fam  (:,:,:,n,nl), &
                                    Time_diag, is, js, 1)
                endif
                if (id_asymdep_fam_column(n,nl)  > 0 ) then
                  used = send_data (id_asymdep_fam_column(n,nl),     &
                               asymdep_fam_col(:,:,n,nl), Time_diag, is, js)
                endif

              end do   ! nl
            end do ! n

            !---- save cmip named fields ----
            do ncmip = 1, size(cmip_family_mapping,1)
              n = cmip_family_mapping(ncmip)
              if (n > 0) then
                if (id_cmipload(ncmip) > 0) then
                  used = send_data (id_cmipload(ncmip), aerosol_fam_col(:,:,n), Time_diag, is, js)
                endif
                if (id_cmipsconc(ncmip) > 0) then
                  used = send_data (id_cmipsconc(ncmip), aerosol_fam(:,:,kerad,n)/deltaz(:,:,kerad), Time_diag, is, js)
                endif
                if (query_cmip_diag_id(ID_cmipconc(ncmip))) then
                  used = send_cmip_data_3d (ID_cmipconc(ncmip), aerosol_fam(:,:,:,n)/deltaz(:,:,:), Time_diag, is, js, 1)
                endif
              endif
            enddo

            if (id_od550aer > 0) then
              used = send_data (id_od550aer, extopdep_fam_col(:,:,naero,nvis), Time_diag, is, js)
            endif
            if (id_abs550aer > 0) then
              used = send_data (id_abs550aer, absopdep_fam_col(:,:,naero,nvis), Time_diag, is, js)
            endif
            if (query_cmip_diag_id(ID_ec550aer)) then
               used = send_cmip_data_3d (ID_ec550aer, extopdep_fam(:,:,:,naero,nvis), Time_diag, is, js, 1)
            endif
            if (id_od550lt1aer > 0) then
              used = send_data (id_od550lt1aer, extopdep_fam_col(:,:,npm25,nvis), Time_diag, is, js)
            endif
            if (id_od870aer > 0) then
              used = send_data (id_od870aer, extopdep_fam_col(:,:,naero,n870), Time_diag, is, js)
            endif
            !----

            deallocate (aerosol_fam)
            deallocate (aerosol_fam_col)
            deallocate (extopdep_fam)
            deallocate (absopdep_fam)
            deallocate (extopdep_fam_col)
            deallocate (absopdep_fam_col)
            if ( Lasymdep ) then
              deallocate (asymdep_fam)
              deallocate (asymdep_fam_col)
              deallocate (sum1)
              deallocate (sum2)
            endif

      endif
    endif
        
!---------------------------------------------------------------------
!    send the user-designated data to diag_manager_mod for processing.
!---------------------------------------------------------------------
        if (id_radswp > 0 ) then
          used = send_data (id_radswp, radswp, Time_diag, is, js, 1)
        endif

        if (id_radp > 0 ) then
          used = send_data (id_radp, radp, Time_diag, is, js, 1)
        endif

        if (id_temp > 0 ) then
          used = send_data (id_temp, temp, Time_diag, is, js, 1)
        endif

        if (id_pressm > 0 ) then
          used = send_data (id_pressm, press, Time_diag, is, js, 1)
        endif

        if (id_phalfm > 0 ) then
          used = send_data (id_phalfm, phalf, Time_diag, is, js, 1)
        endif

        if (id_pfluxm > 0 ) then
          used = send_data (id_pfluxm, pflux, Time_diag, is, js, 1)
        endif

        if (id_rh2o > 0 ) then
          used = send_data (id_rh2o, rh2o, Time_diag, is, js, 1)
        endif

        if (id_cldwater > 0 ) then
          used = send_data (id_cldwater, Cld_spec%cloud_water,  &
                            Time_diag, is, js, 1)
        endif

        if (id_cldice > 0 ) then
          used = send_data (id_cldice, Cld_spec%cloud_ice,  &
                            Time_diag, is, js, 1)
        endif

        if (id_cldarea > 0 ) then
          used = send_data (id_cldarea, Cld_spec%cloud_area,  &
                            Time_diag, is, js, 1)
        endif

        if (id_qo3  > 0 ) then
          used = send_data (id_qo3, qo3, Time_diag, is, js, 1)
        endif

        if (id_qo3v  > 0 ) then
          used = send_data (id_qo3v, 1.0e09*qo3*WTMAIR/WTMOZONE, Time_diag, is, js, 1)
        endif

        !--- save 3D cmip fields
        !--- need log(phalf) for pressure interpolation (what if ptop=0)
        !--- if more fields are added then compute log(phalf) once in the code
        if (query_cmip_diag_id(ID_o3)) then
          used = send_cmip_data_3d (ID_o3, qo3*WTMAIR/WTMOZONE, Time_diag, is, js, 1, phalf=log(phalf))
        endif
        !---

        if (id_qo3_col  > 0 ) then
          used = send_data (id_qo3_col, qo3_col, Time_diag, is, js)
        endif

          if (id_dphalf > 0 ) then
            used = send_data (id_dphalf, dphalf, Time_diag, is, js, 1)
          endif

          if (id_dpflux > 0 ) then
            used = send_data (id_dpflux, dpflux, Time_diag, is, js, 1)
          endif

          if (id_ptop  > 0 ) then
            used = send_data (id_ptop, ptop, Time_diag, is, js)
          endif
        if (Aerosolrad_control%do_aerosol) then
            if (id_sulfate_col_cmip  > 0 ) then
              used = send_data (id_sulfate_col_cmip,      &
                            (96./132.)*aerosol_col(:,:,nso4), Time_diag, is, js)
            endif
            if (id_sulfate_cmip  > 0 ) then
              used = send_data (id_sulfate_cmip,      &
                       (96./132.)*Aerosol%aerosol(:,:,:,nso4)/  &
                                   deltaz(:,:,:), Time_diag, is, js,1)
            endif

            !---- send cmip named fields ----
            if (id_loadso4 > 0) then ! units = vmr
              used = send_data (id_loadso4, (96./132.)*aerosol_col(:,:,nso4), Time_diag, is, js)
            endif
            if (query_cmip_diag_id(ID_concso4)) then
              used = send_cmip_data_3d (ID_concso4, (96./132.)*Aerosol%aerosol(:,:,:,nso4)/deltaz(:,:,:), &
                                        Time_diag, is, js, 1)
            endif
            if (id_sconcso4 > 0) then
              used = send_data (id_sconcso4, (96./132.)*Aerosol%aerosol(:,:,kerad,nso4)/deltaz(:,:,kerad), Time_diag, is, js)
            endif

            if (nsoa > 0) then ! units = mmr
              if (id_loadsoa > 0) then
                used = send_data (id_loadsoa, aerosol_col(:,:,nsoa), Time_diag, is, js)
              endif
              if (query_cmip_diag_id(ID_concsoa)) then
                used = send_cmip_data_3d (ID_concsoa, Aerosol%aerosol(:,:,:,nsoa)/deltaz(:,:,:), Time_diag, is, js, 1)
              endif
              if (id_sconcsoa > 0) then
                used = send_data (id_sconcsoa, Aerosol%aerosol(:,:,kerad,nsoa)/deltaz(:,:,kerad), Time_diag, is, js)
              endif
            endif
            !----

!           if (Aerosolrad_control%do_swaerosol) then
          do n = 1,naerosol
            if (id_aerosol(n)  > 0 ) then
              used = send_data (id_aerosol(n),  &
                                Aerosol%aerosol(:,:,:,n)/deltaz(:,:,:),&
                                Time_diag, is, js, 1)
            endif
            if (id_aerosol_column(n)  > 0 ) then
              used = send_data (id_aerosol_column(n),     &
                                aerosol_col(:,:,n), Time_diag, is, js)
            endif
!           if (Aerosolrad_control%do_swaerosol) then
            do nl=1,N_DIAG_BANDS
              if (id_extopdep(n,nl)  > 0 ) then
                used = send_data (id_extopdep(n,nl),    &
                                  Aerosolrad_diags%extopdep  (:,:,:,n,nl), &
                                  Time_diag, is, js, 1)
              endif
              if (id_extopdep_column(n,nl)  > 0 ) then
                used = send_data (id_extopdep_column(n,nl),     &
                                 extopdep_col(:,:,n,nl), Time_diag, is, js)
              endif
              if (id_absopdep(n,nl)  > 0 ) then
                used = send_data (id_absopdep(n,nl),    &
                                  Aerosolrad_diags%absopdep  (:,:,:,n,nl), &
                                  Time_diag, is, js, 1)
              endif
              if (id_absopdep_column(n,nl)  > 0 ) then
                used = send_data (id_absopdep_column(n,nl),     &
                                 absopdep_col(:,:,n,nl), Time_diag, is, js)
              endif
              if (id_asymdep(n,nl)  > 0 ) then
                used = send_data (id_asymdep(n,nl),    &
                                  Aerosolrad_diags%asymdep  (:,:,:,n,nl), &
                                  Time_diag, is, js, 1)
              endif
              if (id_asymdep_column(n,nl)  > 0 ) then
                used = send_data (id_asymdep_column(n,nl),     &
                                 asymdep_col(:,:,n,nl), Time_diag, is, js)
              endif
!           endif
            end do
          end do
          if (Aerosolrad_control%volcanic_lw_aerosols) then
!           co_indx = size(Aerosolrad_diags%lw_ext,4)
            co_indx = 5
            if (id_lwext_vlcno(1)  > 0 ) then
              used = send_data (id_lwext_vlcno(1),     &
                         Aerosolrad_diags%lw_ext(:,:,:,co_indx)*deltaz(:,:,:), &
                         Time_diag, is, js,1)
            endif
            if (id_lw_xcoeff_vlcno(1)  > 0 ) then
              used = send_data (id_lw_xcoeff_vlcno(1),     &
                            Aerosolrad_diags%lw_ext(:,:,:,co_indx),  &
                            Time_diag, is, js,1)
            endif
            if (id_lwssa_vlcno(1)  > 0 ) then
              used = send_data (id_lwssa_vlcno(1),     &
                            Aerosolrad_diags%lw_ssa(:,:,:,co_indx),  &
                            Time_diag, is, js,1)
            endif
            if (id_lwasy_vlcno(1)  > 0 ) then
              used = send_data (id_lwasy_vlcno(1),     &
                           Aerosolrad_diags%lw_asy(:,:,:,co_indx),  &
                           Time_diag, is, js,1)
            endif
!           bnd_indx = 4
            bnd_indx = 6
            if (id_lwext_vlcno(2)  > 0 ) then
              used = send_data (id_lwext_vlcno(2),     &
                               Aerosolrad_diags%lw_ext(:,:,:,bnd_indx)*deltaz(:,:,:), &
                               Time_diag, is, js,1)
            endif
            if (id_lw_xcoeff_vlcno(2)  > 0 ) then
              used = send_data (id_lw_xcoeff_vlcno(2),     &
                               Aerosolrad_diags%lw_ext(:,:,:,bnd_indx),  &
                               Time_diag, is, js,1)
            endif
            if (id_lwssa_vlcno(2)  > 0 ) then
              used = send_data (id_lwssa_vlcno(2),     &
                               Aerosolrad_diags%lw_ssa(:,:,:,bnd_indx), &
                               Time_diag, is, js,1)
            endif
            if (id_lwasy_vlcno(2)  > 0 ) then
              used = send_data (id_lwasy_vlcno(2),     &
                              Aerosolrad_diags%lw_asy(:,:,:,bnd_indx), &
                              Time_diag, is, js,1)
            endif
            do nv=1,2
              if (id_lw_extopdep_vlcno_column(nv)  > 0 ) then
                used = send_data (id_lw_extopdep_vlcno_column(nv),     &
                                  lw_extopdep_vlcno_col(:,:,nv),  &
                                  Time_diag, is, js)
              endif
              if (id_lw_absopdep_vlcno_column(nv)  > 0 ) then
                used = send_data (id_lw_absopdep_vlcno_column(nv),     &
                                  lw_absopdep_vlcno_col(:,:,nv),  &
                                  Time_diag, is, js)
              endif
            end do
            deallocate (lw_absopdep_vlcno_col, lw_extopdep_vlcno_col)
          endif
          if (Aerosolrad_control%volcanic_sw_aerosols) then
            if (id_swext_vlcno(1)  > 0 ) then
              used = send_data (id_swext_vlcno(1),     &
                                Aerosolrad_diags%sw_ext(:,:,:,w550_band_indx)*deltaz(:,:,:), &
                                Time_diag, is, js,1)
            endif
            if (id_sw_xcoeff_vlcno(1)  > 0 ) then
              used = send_data (id_sw_xcoeff_vlcno(1),     &
                               Aerosolrad_diags%sw_ext(:,:,:,w550_band_indx),  &
                               Time_diag, is, js,1)
            endif
            if (id_swssa_vlcno(1)  > 0 ) then
              used = send_data (id_swssa_vlcno(1),     &
                                 Aerosolrad_diags%sw_ssa(:,:,:,w550_band_indx), &
                                 Time_diag, is, js,1)
            endif
            if (id_swasy_vlcno(1)  > 0 ) then
              used = send_data (id_swasy_vlcno(1),     &
                                Aerosolrad_diags%sw_asy(:,:,:,w550_band_indx), &
                                Time_diag, is, js,1)
            endif
            if (id_swheat_vlcno  > 0 ) then
              v_heat = Aerosolrad_diags%sw_heating_vlcno(:,:,:)/float(1) ! for repro
              used = send_data (id_swheat_vlcno ,    &
!                               Aerosolrad_diags%sw_heating_vlcno(:,:,:), &
                                v_heat(:,:,:), &
                                Time_diag, is, js, 1)
            endif
            if (id_swext_vlcno(2)  > 0 ) then
              used = send_data (id_swext_vlcno(2),     &
                               Aerosolrad_diags%sw_ext(:,:,:,one_micron_indx)*deltaz(:,:,:), &
                               Time_diag, is, js,1)
            endif
            if (id_sw_xcoeff_vlcno(2)  > 0 ) then
              used = send_data (id_sw_xcoeff_vlcno(2),     &
                               Aerosolrad_diags%sw_ext(:,:,:,one_micron_indx),  &
                               Time_diag, is, js,1)
            endif
            if (id_swssa_vlcno(2)  > 0 ) then
              used = send_data (id_swssa_vlcno(2),     &
                                Aerosolrad_diags%sw_ssa(:,:,:,one_micron_indx),  &
                                Time_diag, is, js,1)
            endif
            if (id_swasy_vlcno(2)  > 0 ) then
              used = send_data (id_swasy_vlcno(2),     &
                                Aerosolrad_diags%sw_asy(:,:,:,one_micron_indx), &
                                Time_diag, is, js,1)
            endif
            do nv=1,2
              if (id_extopdep_vlcno_column(nv)  > 0 ) then
                used = send_data (id_extopdep_vlcno_column(nv),     &
                                  extopdep_vlcno_col(:,:,nv),  &
                                  Time_diag, is, js)
              endif
              if (id_absopdep_vlcno_column(nv)  > 0 ) then
                used = send_data (id_absopdep_vlcno_column(nv),     &
                                  absopdep_vlcno_col(:,:,nv),  &
                                  Time_diag, is, js)
              endif
            end do
            deallocate (absopdep_vlcno_col, extopdep_vlcno_col)
          endif
          deallocate (aerosol_col)
!         if (Aerosolrad_control%do_swaerosol) then
            deallocate (extopdep_col)
            deallocate (absopdep_col)
            if (allocated (asymdep_col)) deallocate (asymdep_col)
!         endif
        endif

        if (id_cmxolw > 0 ) then
          used = send_data (id_cmxolw, cmxolw, Time_diag, is, js, 1)
        endif

        if (id_crndlw > 0 ) then
          used = send_data (id_crndlw, crndlw, Time_diag, is, js, 1)
        endif

        if (id_flxnet > 0 ) then
          used = send_data (id_flxnet, flxnet, Time_diag, is, js, 1)
        endif

        if (id_fsw    > 0 ) then
          used = send_data (id_fsw   , fsw   , Time_diag, is, js, 1)
        endif

        if (id_ufsw   > 0 ) then
          used = send_data (id_ufsw  , ufsw  , Time_diag, is, js, 1)
        endif

        if (id_dfsw   > 0 ) then
          used = send_data (id_dfsw  , ufsw-fsw  , Time_diag, is, js, 1)
        endif

        if (id_psj    > 0 ) then
          used = send_data (id_psj   , psj   , Time_diag, is, js)
        endif

        if (id_tmpsfc > 0 ) then
          used = send_data (id_tmpsfc, tmpsfc, Time_diag, is, js)
        endif

        if (id_cvisrfgd_dir > 0 ) then
          used = send_data (id_cvisrfgd_dir, asfc_vis_dir, Time_diag, is, js)
        endif
 
        if (id_cvisrfgd_dif > 0 ) then
          used = send_data (id_cvisrfgd_dif, asfc_vis_dif, Time_diag, is, js)
        endif

        if (id_cirrfgd_dir > 0 ) then
          used = send_data (id_cirrfgd_dir , asfc_nir_dir, Time_diag, is,js)
        endif
 
        if (id_cirrfgd_dif > 0 ) then
          used = send_data (id_cirrfgd_dif , asfc_nir_dif, Time_diag, is,js)
        endif

     do n=1, 4
        if (id_sw_bdyflx(n) > 0 ) then
          used = send_data (id_sw_bdyflx(n) , bdy_flx_mean(:,:,n),&
                            Time_diag, is,js)
        endif
     end do
! change by Xianglei Huang, output more lw_bdyflx
     do n=1, 7
        if (id_lw_bdyflx(n) > 0 ) then
          used = send_data (id_lw_bdyflx(n) , Lw_output%bdy_flx(:,:,n),&
                            Time_diag, is,js)
        endif
     end do
     if (id_bT_Aqua31 > 0) then
        !Assumes band 3 is 800--900
        !Aqua31 recieves between 885 and 925 cm**-1, using band 3
        !Computation from http://pds-atmospheres.nmsu.edu/education_and_outreach/encyclopedia/planck_function.htm
        !Note their planck function is in mW and bdy_flx is is W
        !850 cm**-1 gives a characteristic wavenumber for this band, which is 100 cm**-1 wide
        !Also, we divide by pi to convert the omni-directional flux into intensity
        bT = 0.0000119*(850.**3)/((Lw_output%bdy_flx(:,:,3)/100.)*1000./pi + 1.)
        bT = 1.438833*850./log(bT)
        used = send_data ( id_bT_Aqua31, bT, Time_diag, is, js ) 
     endif

        if (Rad_control%do_totcld_forcing) then
     do n=1, 4
        if (id_sw_bdyflx_clr(n) > 0 ) then
          used = send_data (id_sw_bdyflx_clr(n) ,   &
                            bdy_flx_clr_mean(:,:,n),&
                            Time_diag, is,js)
        endif
     end do
! change by Xianglei Huang, output more lw_bdyflx_clr
     do n=1, 7
        if (id_lw_bdyflx_clr(n) > 0 ) then
          used = send_data (id_lw_bdyflx_clr(n) ,   &
                            Lw_output%bdy_flx_clr(:,:,n),&
                            Time_diag, is,js)
        endif
     end do

          if (id_radswpcf > 0 ) then
            used = send_data (id_radswpcf, radswpcf, Time_diag,   &
                              is, js, 1)
          endif

          if (id_radpcf > 0 ) then
            used = send_data (id_radpcf, radpcf, Time_diag, is, js, 1)
          endif

          if (id_flxnetcf > 0 ) then
            used = send_data (id_flxnetcf, flxnetcf, Time_diag,   &
                              is, js, 1)
          endif

          if (id_fswcf  > 0 ) then
            used = send_data (id_fswcf , fswcf , Time_diag, is, js, 1)
          endif

          if (id_ufswcf  > 0 ) then
            used = send_data (id_ufswcf , ufswcf , Time_diag, is, js, 1)
          endif

          if (id_dfswcf  > 0 ) then
            used = send_data (id_dfswcf , ufswcf-fswcf , Time_diag, is, js, 1)
          endif
        endif

      endif    ! (Time_diag > Time)
      endif    ! (write_data_file)

!------------------------------------------------------------------


end subroutine write_rad_output_file



!#####################################################################
! <SUBROUTINE NAME="rad_output_file_end">
!  <OVERVIEW>
!   rad_output_file_end is the destructor for rad_output_file_mod
!  </OVERVIEW>
!  <DESCRIPTION>
!   rad_output_file_end is the destructor for rad_output_file_mod
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rad_output_file_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine rad_output_file_end

!-------------------------------------------------------------------
!    rad_output_file_end is the destructor for rad_output_file_mod.
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_output_file_mod', &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized= .false. 

!----------------------------------------------------------------------

end subroutine rad_output_file_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#################################################################
! <SUBROUTINE NAME="register_fields">
!  <OVERVIEW>
!   register_fields send the relevant information concerning the 
!    user-desired output fields to diag_manager_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   register_fields send the relevant information concerning the 
!    user-desired output fields to diag_manager_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call register_fields (Time, axes, nfilds, names)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="nfields" TYPE="integer">
!   number of aerosol fields
!  </IN>
!  <IN NAME="names" TYPE="character">
!   names of aerosol fields
!  </IN>
! </SUBROUTINE>
!
subroutine register_fields (Time, axes, nfields, names, family_names, &
                            do_totcld_forcing, volcanic_lw_aerosols, &
                            volcanic_sw_aerosols)

!--------------------------------------------------------------------
!    register_fields send the relevant information concerning the 
!    user-desired output fields to diag_manager_mod.
!--------------------------------------------------------------------

type(time_type),                intent(in) :: Time
integer, dimension(4),          intent(in) :: axes
integer,                        intent(in) :: nfields
character(len=*), dimension(:), intent(in) :: names, family_names
logical,                        intent(in) :: do_totcld_forcing
logical,                        intent(in) :: volcanic_lw_aerosols
logical,                        intent(in) :: volcanic_sw_aerosols

!--------------------------------------------------------------------
!  intent(in) variables:
!
!       Time      current time [ time_type(days, seconds) ]
!       axes      diagnostic variable axes for netcdf files
!       nfields   number of active aerosol species
!       names     names of active aerosol species
!       family_names  
!                 names of active aerosol families
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      character(len=64), dimension(:), allocatable ::   & 
                                               aerosol_column_names, &
                                              extopdep_column_names, &
                                              absopdep_column_names, &
                                                 extopdep_names, &
                                                 absopdep_names, &
                                             asymdep_column_names, &
                                             asymdep_names, &
                                             asymdep_fam_column_names, &
                                             asymdep_fam_names, &
                                          aerosol_fam_column_names, &
                                         extopdep_fam_column_names, &
                                        absopdep_fam_column_names, &
                                        extopdep_fam_names, &
                                        absopdep_fam_names
      integer, dimension(4)    :: bxes
      integer                  :: n, nl
      integer                  :: nfamilies
      real                     :: trange(2)
      integer                  :: ncmip

!---------------------------------------------------------------------
!   local variables:
!
!       aerosol_column_names
!       bxes                   diagnostic variable axes with elements 
!                              (1:3) valid for variables defined at
!                              flux levels
!       n
!
!--------------------------------------------------------------------
 
!-------------------------------------------------------------------
!    define variable axis array with elements (1:3) valid for variables
!    defined at flux levels.
!-------------------------------------------------------------------
      bxes(1:2) = axes(1:2)
      bxes(3) = axes(4)
      bxes(4) = axes(4)
      trange =(/ 100., 400. /)

!---------------------------------------------------------------------
!    register the potential diagnostic variables from this module.
!    the variables will actually be saved and then output only if
!    they are activated through the input diag_file.
!---------------------------------------------------------------------
      id_radswp = &
         register_diag_field (mod_name, 'radswp', axes(1:3), Time, &
                          'temperature tendency for SW radiation', &
                          'deg_K/sec', missing_value=missing_value)

      id_radp = &
         register_diag_field (mod_name, 'radp', axes(1:3), Time, &
                          'temperature tendency for radiation', &
                          'deg_K/sec', missing_value=missing_value)

      id_temp   = &
         register_diag_field (mod_name, 'temp', axes(1:3), Time, &
                          'temperature field', &
                          'deg_K', missing_value=trange(1), &
                          range=trange)

      id_pressm  = &
         register_diag_field (mod_name, 'pressm', axes(1:3), Time, &
                           'model level pressure', &
                           'Pa', missing_value=missing_value)
 
      id_phalfm  = &
         register_diag_field (mod_name, 'phalfm', bxes(1:3), Time, &
                           'model interface level pressure', &
                           'Pa', missing_value=missing_value)

      id_pfluxm  = &
          register_diag_field (mod_name, 'pfluxm', bxes(1:3), Time, &
                           'radiation flux level pressures', &
                           'Pa', missing_value=missing_value)

        id_dpflux  = &
          register_diag_field (mod_name, 'dpflux', axes(1:3), Time,  &
                           'radiation flux layer thickness', &
                           'hPa', missing_value=missing_value)
   
        id_dphalf  = &
          register_diag_field (mod_name, 'dphalf', axes(1:3), Time,  &
                           'radiation model layer thickness', &
                           'hPa', missing_value=missing_value)
   
        id_ptop  = &
          register_diag_field (mod_name, 'ptop', axes(1:2), Time, &
                           'pressure at model top', &
                           'hPa', missing_value=missing_value)
 

      id_rh2o   = &
         register_diag_field (mod_name, 'rh2o', axes(1:3), Time, &
                          'water vapor mixing ratio', &
                          'kg/kg', missing_value=missing_value)

      id_cldwater = &
         register_diag_field (mod_name, 'cloud_water', axes(1:3), Time,&
                          'cloud water specific humidity', &
                          'kg/kg', missing_value=missing_value)

      id_cldice = &
         register_diag_field (mod_name, 'cloud_ice', axes(1:3), Time,&
                          'cloud ice specific humidity', &
                          'kg/kg', missing_value=missing_value)

      id_cldarea = &
         register_diag_field (mod_name, 'cloud_area', axes(1:3), Time,&
                          'cloud fractional area', &
                          'fraction', missing_value=missing_value)

      id_qo3    = &
         register_diag_field (mod_name, 'qo3', axes(1:3), Time, &
                          'ozone mixing ratio', &
                          'kg/kg', missing_value=missing_value)

      id_qo3v    = &
         register_diag_field (mod_name, 'qo3v', axes(1:3), Time, &
                          'ozone mole fraction', &
                          '1.e-9', missing_value=missing_value)

      ID_o3    = register_cmip_diag_field_3d (mod_name, 'o3', Time, &
                          'Ozone Volume Mixing Ratio', 'mol mol-1', &
                        standard_name = 'mole_fraction_of_ozone_in_air')
              
      id_qo3_col = &
         register_diag_field (mod_name, 'qo3_col', axes(1:2), Time, &
                          'ozone column', &
                          'DU', missing_value=missing_value)

!--------------------------------------------------------------------
!    allocate space for and save aerosol name information.
!--------------------------------------------------------------------
      if (nfields /= 0) then
        naerosol = nfields
        allocate (id_aerosol(naerosol))
        allocate (id_aerosol_column(naerosol)) 
        allocate (aerosol_column_names(naerosol))
        id_sulfate_col_cmip = &
             register_diag_field (mod_name, 'sulfate_col_cmip',  &
                            axes(1:2), Time, 'sulfate_col_cmip',&
                                  'kg/m2', missing_value=missing_value)
        id_sulfate_cmip = &
             register_diag_field (mod_name, 'sulfate_cmip',  &
                            axes(1:3), Time, 'sulfate_cmip',&
                                  'kg/m3', missing_value=missing_value)

        nso4 = 0
        nsoa = 0
        do n = 1, naerosol                           
          aerosol_column_names(n) = TRIM(names(n) ) // "_col"
          if (lowercase(TRIM(names(n))) == 'so4')      nso4 = n
          if (lowercase(TRIM(names(n))) == 'omphilic') nsoa = n
          ! DEBUG
         !call error_mesg('rad_output_file_mod','tracer_name='//lowercase(TRIM(names(n)))//', nsoa='//TRIM(string(nsoa)),NOTE)
        end do

        ! register cmip named fields
        if (nso4 > 0) then
          id_loadso4 = register_cmip_diag_field_2d (mod_name, 'loadso4',  &
                                         Time, 'Load of SO4', 'kg m-2', &
                        standard_name='atmosphere_mass_content_of_sulfate_dry_aerosol')
          ID_concso4 = register_cmip_diag_field_3d (mod_name, 'concso4', Time, &
                            'Concentration of SO4', 'kg m-3', &
                        standard_name='mass_concentration_of_sulfate_dry_aerosol_in_air')
          id_sconcso4 = register_cmip_diag_field_2d (mod_name, 'sconcso4', Time, &
                                     'Surface Concentration of SO4', 'kg m-3', &
                        standard_name='mass_concentration_of_sulfate_dry_aerosol_in_air')
        else
          id_loadso4 = 0
          id_sconcso4 = 0
        endif

        if (nsoa > 0) then
          id_loadsoa = register_cmip_diag_field_2d (mod_name, 'loadsoa', Time,  &
                                    'Load of Dry Aerosol Secondary Organic Matter', 'kg m-2', &
                    standard_name='atmosphere_mass_content_of_secondary_particulate_organic_matter_dry_aerosol')
          id_sconcsoa = register_cmip_diag_field_2d (mod_name, 'sconcsoa', Time,  &
                               'Surface Concentration of Dry Aerosol Secondary Organic Matter', 'kg m-3', &
                          standard_name='mass_concentration_of_secondary_particulate_organic_matter_dry_aerosol_in_air')
          ID_concsoa  = register_cmip_diag_field_3d (mod_name, 'concsoa', Time,  &
                               'Concentration of Dry Aerosol Secondary Organic Matter', 'kg m-3', &
                          standard_name='mass_concentration_of_secondary_particulate_organic_matter_dry_aerosol_in_air')
        else
          id_loadsoa = 0
          id_sconcsoa = 0
        endif

        !---- register tracer/aerosol diagnostics ----
        do n = 1,naerosol
          id_aerosol(n)    = &
             register_diag_field (mod_name, TRIM(names(n)), axes(1:3), &
                                  Time, TRIM(names(n)),&
                                  'kg/m3', missing_value=missing_value)
          id_aerosol_column(n)    = &
             register_diag_field (mod_name,   &
                      TRIM(aerosol_column_names(n)), axes(1:2), Time, &
                      TRIM(aerosol_column_names(n)), &
                      'kg/m2', missing_value=missing_value)
        end do
        deallocate (aerosol_column_names)

        allocate (extopdep_names(naerosol))
        allocate (extopdep_column_names(naerosol))
        allocate (absopdep_names(naerosol))
        allocate (absopdep_column_names(naerosol))
        allocate (id_extopdep(naerosol, N_DIAG_BANDS))
        allocate (id_extopdep_column(naerosol, N_DIAG_BANDS))
        allocate (id_absopdep(naerosol, N_DIAG_BANDS))
        allocate (id_absopdep_column(naerosol, N_DIAG_BANDS))
        allocate (id_asymdep(naerosol, N_DIAG_BANDS))
        allocate (id_asymdep_column(naerosol, N_DIAG_BANDS))
        allocate (asymdep_names(naerosol))
        allocate (asymdep_column_names(naerosol))
     do nl=1,N_DIAG_BANDS
        do n = 1,naerosol                           
          extopdep_names(n) =   &
                TRIM(names(n) ) // "_exopdep" // TRIM(band_suffix(nl))
          extopdep_column_names(n) =   &
             TRIM(names(n) ) // "_exopdep_col" // TRIM(band_suffix(nl))
          absopdep_names(n) =   &
             TRIM(names(n) ) // "_abopdep" // TRIM(band_suffix(nl))
          absopdep_column_names(n) =   &
             TRIM(names(n) ) // "_abopdep_col" // TRIM(band_suffix(nl))
          asymdep_names(n) =   &
             TRIM(names(n) ) // "_asymdep" // TRIM(band_suffix(nl))
          asymdep_column_names(n) =   &
             TRIM(names(n) ) // "_asymdep_col" // TRIM(band_suffix(nl))
        end do
        do n = 1,naerosol
          id_extopdep(n,nl)    = &
             register_diag_field (mod_name, TRIM(extopdep_names(n)), axes(1:3), &
                                  Time, TRIM(extopdep_names(n)),&
                                  'dimensionless', missing_value=missing_value)
          id_extopdep_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(extopdep_column_names(n)), axes(1:2), Time, &
                      TRIM(extopdep_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
          id_absopdep(n,nl)    = &
             register_diag_field (mod_name, TRIM(absopdep_names(n)), axes(1:3), &
                                  Time, TRIM(absopdep_names(n)),&
                                  'dimensionless', missing_value=missing_value)
          id_absopdep_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(absopdep_column_names(n)), axes(1:2), Time, &
                      TRIM(absopdep_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
          id_asymdep(n,nl)    = &
             register_diag_field (mod_name, TRIM(asymdep_names(n)), axes(1:3),&
                                  Time, TRIM(asymdep_names(n)),&
                                  'dimensionless', missing_value=missing_value)
          id_asymdep_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(asymdep_column_names(n)), axes(1:2), Time, &
                      TRIM(asymdep_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
        end do
      end do
        deallocate (extopdep_names)
        deallocate (extopdep_column_names)
        deallocate (absopdep_names)
        deallocate (absopdep_column_names)
        deallocate (asymdep_names)
        deallocate (asymdep_column_names)
      endif

      if (size(family_names(:)) /= 0) then
        nfamilies = size(family_names(:))
        allocate (id_aerosol_fam(nfamilies))
        allocate (id_aerosol_fam_column(nfamilies)) 
        allocate (aerosol_fam_column_names(nfamilies))
        do n=1,nfamilies      
          aerosol_fam_column_names(n) = TRIM(family_names(n) ) // "_col"
        end do
        cmip_family_mapping = 0 ! mapping between family index and cmip diag fields
        do n = 1,nfamilies
          id_aerosol_fam(n)    = &
             register_diag_field (mod_name, TRIM(family_names(n)), axes(1:3), &
                                  Time, TRIM(family_names(n)),&
                                  'kg/m3', missing_value=missing_value)
          id_aerosol_fam_column(n)    = &
             register_diag_field (mod_name,   &
                      TRIM(aerosol_fam_column_names(n)), axes(1:2), Time, &
                      TRIM(aerosol_fam_column_names(n)), &
                      'kg/m2', missing_value=missing_value)

          !---- cmip diagnostic fields (load, sconc, conc) ----
          ncmip = 0
          if (TRIM(family_names(n)) .eq. 'organic_carbon') ncmip = 1
          if (TRIM(family_names(n)) .eq. 'POA')            ncmip = 2
          if (TRIM(family_names(n)) .eq. 'SOA' .and. &
                               nsoa .eq. 0 )               ncmip = 3
          if (TRIM(family_names(n)) .eq. 'black_carbon')   ncmip = 4
          if (TRIM(family_names(n)) .eq. 'dust')           ncmip = 5
          if (TRIM(family_names(n)) .eq. 'sea_salt')       ncmip = 6
          call error_mesg('rad_output_file_mod','family_name='//TRIM(family_names(n))//', ncmip='//TRIM(string(ncmip)),NOTE)
          if (ncmip > 0) then
            cmip_family_mapping(ncmip) = n
            id_cmipload(ncmip) = register_cmip_diag_field_2d (mod_name, 'load'//TRIM(cmip_names(ncmip)), &
                                                Time, 'Load of '//TRIM(cmip_longnames(ncmip)), 'kg m-2', &
                          standard_name='atmosphere_mass_content_of_'//TRIM(cmip_stdnames(ncmip))//'_dry_aerosol')
            id_cmipsconc(ncmip) = register_cmip_diag_field_2d (mod_name, 'sconc'//TRIM(cmip_names(ncmip)), &
                                 Time, 'Surface Concentration of '//TRIM(cmip_longnames(ncmip)), 'kg m-3', &
                          standard_name='mass_concentration_of_'//TRIM(cmip_stdnames(ncmip))//'_dry_aerosol_in_air')
            ID_cmipconc(ncmip) = register_cmip_diag_field_3d (mod_name, 'conc'//TRIM(cmip_names(ncmip)), &
                                 Time, 'Concentration of '//TRIM(cmip_longnames(ncmip)), 'kg m-3', &
                          standard_name='mass_concentration_of_'//TRIM(cmip_stdnames(ncmip))//'_dry_aerosol_in_air')
          endif
          !----

        end do
        deallocate (aerosol_fam_column_names)

        allocate (id_extopdep_fam(nfamilies, N_DIAG_BANDS))
        allocate (id_extopdep_fam_column(nfamilies, N_DIAG_BANDS))
        allocate (id_absopdep_fam(nfamilies, N_DIAG_BANDS))
        allocate (id_absopdep_fam_column(nfamilies, N_DIAG_BANDS))
        allocate (extopdep_fam_names(nfamilies))
        allocate (extopdep_fam_column_names(nfamilies))
        allocate (absopdep_fam_names(nfamilies))
        allocate (absopdep_fam_column_names(nfamilies))
        allocate (id_asymdep_fam(nfamilies, N_DIAG_BANDS))
        allocate (id_asymdep_fam_column(nfamilies, N_DIAG_BANDS))
        allocate (asymdep_fam_names(nfamilies))
        allocate (asymdep_fam_column_names(nfamilies))

        ! indices to aerosol, pm2.5, vis, and 870 bands
        naero = 0; npm25 = 0; nvis = 0; n870 = 0
        do n = 1, nfamilies
          if (TRIM(family_names(n)) .eq. 'aerosol') naero = n
          if (TRIM(family_names(n)) .eq. 'pm2.5')   npm25 = n
        enddo
        do nl = 1, N_DIAG_BANDS
          if (TRIM(band_suffix(nl)) .eq. '_vis') nvis = nl
          if (TRIM(band_suffix(nl)) .eq. '_870') n870 = nl
        enddo

   do nl=1,N_DIAG_BANDS
        do n=1,nfamilies      
          extopdep_fam_names(n) =   &
           TRIM(family_names(n) ) // "_exopdep" // TRIM(band_suffix(nl))
          extopdep_fam_column_names(n) =   &
       TRIM(family_names(n) ) // "_exopdep_col" // TRIM(band_suffix(nl))
          absopdep_fam_names(n) =   &
          TRIM(family_names(n) ) // "_abopdep" // TRIM(band_suffix(nl))
          absopdep_fam_column_names(n) =  &
       TRIM(family_names(n) ) // "_abopdep_col" // TRIM(band_suffix(nl))
          asymdep_fam_names(n) =   &
          TRIM(family_names(n) ) // "_asymdep" // TRIM(band_suffix(nl))
          asymdep_fam_column_names(n) =  &
       TRIM(family_names(n) ) // "_asymdep_col" // TRIM(band_suffix(nl))
        end do
        do n = 1,nfamilies
          id_extopdep_fam(n,nl)    = &
             register_diag_field (mod_name, TRIM(extopdep_fam_names(n)), axes(1:3), &
                                  Time, TRIM(extopdep_fam_names(n)),&
                                  'dimensionless', missing_value=missing_value)
          id_extopdep_fam_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(extopdep_fam_column_names(n)), axes(1:2), Time, &
                      TRIM(extopdep_fam_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
          id_absopdep_fam(n,nl)    = &
             register_diag_field (mod_name, TRIM(absopdep_fam_names(n)), axes(1:3), &
                                  Time, TRIM(absopdep_fam_names(n)),&
                                  'dimensionless', missing_value=missing_value)
          id_absopdep_fam_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(absopdep_fam_column_names(n)), axes(1:2), Time, &
                      TRIM(absopdep_fam_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
          id_asymdep_fam(n,nl)    = &
             register_diag_field(mod_name,TRIM(asymdep_fam_names(n)),axes(1:3),&
                      Time, TRIM(asymdep_fam_names(n)),&
                     'dimensionless', missing_value=missing_value)
          id_asymdep_fam_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(asymdep_fam_column_names(n)), axes(1:2), Time, &
                      TRIM(asymdep_fam_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
        end do
   end do

        !---- register cmip fields ----
        id_od550aer = 0; id_abs550aer = 0; id_od550lt1aer = 0; id_od870aer = 0
        if (naero > 0 .and. nvis > 0) then
          id_od550aer = register_cmip_diag_field_2d (mod_name, 'od550aer', Time, &
                            'Ambient Aerosol Optical Thickness at 550 nm', '1.0', &
                            standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol_particles')
          id_abs550aer = register_cmip_diag_field_2d (mod_name, 'abs550aer', Time, &
                            'Ambient Aerosol Absorption Optical Thickness at 550 nm', '1.0', &
                            standard_name='atmosphere_absorption_optical_thickness_due_to_ambient_aerosol')
          ID_ec550aer = register_cmip_diag_field_3d (mod_name, 'ec550aer', Time, &
                            'Ambient Aerosol Extinction at 550 nm', 'm-1', &
                            standard_name='volume_extinction_coefficient_in_air_due_to_ambient_aerosol')
        endif
        if (npm25 > 0 .and. nvis > 0) then
          id_od550lt1aer = register_cmip_diag_field_2d (mod_name, 'od550lt1aer', Time, &
                            'Ambient Fine Aerosol Optical Thickness at 550 nm', '1.0', &
                            standard_name='atmosphere_optical_thickness_due_to_pm1_ambient_aerosol')
        endif
        if (naero > 0 .and. n870 > 0) then
          id_od870aer = register_cmip_diag_field_2d (mod_name, 'od870aer', Time, &
                            'Ambient Aerosol Optical Thickness at 870 nm', '1.0', &
                            standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol_particles')
        endif
        !----

        deallocate (extopdep_fam_names)
        deallocate (extopdep_fam_column_names)
        deallocate (absopdep_fam_names)
        deallocate (absopdep_fam_column_names)
        deallocate (asymdep_fam_names)
        deallocate (asymdep_fam_column_names)

      endif
      
        if (volcanic_lw_aerosols) then
          id_lw_extopdep_vlcno_column(1) = &
             register_diag_field (mod_name,   &
                    'lw_b5_extopdep_vlcno_c', axes(1:2), Time, &
                    'lw 900-990 band column volcanic extopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_absopdep_vlcno_column(1)    = &
             register_diag_field (mod_name,   &
                    'lw_b5_absopdep_vlcno_c', axes(1:2), Time, &
                    'lw 900-990 band column volcanic absopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_lwext_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'bnd5_extopdep_vlcno', axes(1:3), Time, &
                    '900-990 band volcanic lw extopdep  ',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_xcoeff_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'bnd5_lwext_vlcno', axes(1:3), Time, &
                    '900-990 band volcanic lw extinction',  &
                      'meter**(-1)  ', missing_value=missing_value)
          id_lwssa_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'bnd5_lwssa_vlcno', axes(1:3), Time, &
                    '900-990   band volcanic lw scattering albedo', &
                      'dimensionless', missing_value=missing_value)
          id_lwasy_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'bnd5_lwasy_vlcno', axes(1:3), Time, &
                    '900-990 band volcanic lw asymmetry',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_extopdep_vlcno_column(2) = &
             register_diag_field (mod_name,   &
                    'bnd6_extopdep_vlcno_c', axes(1:2), Time, &
                    '990-1070 column volcanic extopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_absopdep_vlcno_column(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_absopdep_vlcno_c', axes(1:2), Time, &
                    '990-1070 column volcanic absopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_lwext_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_extopdep_vlcno', axes(1:3), Time, &
                    '990-1070 volcanic lw extopdep  ',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_xcoeff_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_lwext_vlcno', axes(1:3), Time, &
                    '990-1070 volcanic lw extinction',  &
                      'meter**(-1)  ', missing_value=missing_value)
          id_lwssa_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_lwssa_vlcno', axes(1:3), Time, &
                    '990-1070 volcanic lw scattering albedo',  &
                      'dimensionless', missing_value=missing_value)
          id_lwasy_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_lwasy_vlcno', axes(1:3), Time, &
                    '990-1070 volcanic lw asymmetry',  &
                      'dimensionless', missing_value=missing_value)
        endif

        if (volcanic_sw_aerosols) then
          id_extopdep_vlcno_column(1) = &
             register_diag_field (mod_name,   &
                    'vis_extopdep_vlcno_c', axes(1:2), Time, &
                    'visband column volcanic extopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_absopdep_vlcno_column(1)    = &
             register_diag_field (mod_name,   &
                    'vis_absopdep_vlcno_c', axes(1:2), Time, &
                    'visband column volcanic absopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_swext_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'visband_swextopdep_vlcno', axes(1:3), Time, &
                    'visband volcanic sw extopdep  ',  &
                      'dimensionless', missing_value=missing_value)
          id_sw_xcoeff_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'visband_swext_vlcno', axes(1:3), Time, &
                    'visband volcanic sw extinction',  &
                      'meter**(-1)', missing_value=missing_value)
          id_swssa_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'visband_swssa_vlcno', axes(1:3), Time, &
                    'visband volcanic sw scattering albedo',  &
                      'dimensionless', missing_value=missing_value)
          id_swasy_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'visband_swasy_vlcno', axes(1:3), Time, &
                    'visband volcanic sw asymmetry',  &
                      'dimensionless', missing_value=missing_value)
          id_extopdep_vlcno_column(2) = &
             register_diag_field (mod_name,   &
                    'nir_extopdep_vlcno_c', axes(1:2), Time, &
                    'nirband column volcanic extopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_absopdep_vlcno_column(2)    = &
             register_diag_field (mod_name,   &
                    'nir_absopdep_vlcno_c', axes(1:2), Time, &
                    'nirband column volcanic absopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_swext_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'nirband_swextopdep_vlcno', axes(1:3), Time, &
                    'nirband volcanic sw extopdep  ',  &
                      'dimensionless', missing_value=missing_value)
          id_sw_xcoeff_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'nirband_swext_vlcno', axes(1:3), Time, &
                    'nirband volcanic sw extinction',  &
                      'meter**(-1)', missing_value=missing_value)
          id_swssa_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'nirband_swssa_vlcno', axes(1:3), Time, &
                    'nirband volcanic sw scattering albedo',  &
                      'dimensionless', missing_value=missing_value)
          id_swasy_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'nirband_swasy_vlcno', axes(1:3), Time, &
                    'nirband volcanic sw asymmetry',  &
                      'dimensionless', missing_value=missing_value)
          id_swheat_vlcno    = &
             register_diag_field (mod_name,   &
                    'sw_heating_vlcno', axes(1:3), Time, &
                    'sw heating due to vlcnic aero',  &
                      'deg K per day', missing_value=missing_value)
        endif


      id_cmxolw = &
         register_diag_field (mod_name, 'cmxolw', axes(1:3), Time, &
                          'maximum overlap cloud amount', &
                          'percent', missing_value=missing_value)

      id_crndlw = &
         register_diag_field (mod_name, 'crndlw', axes(1:3), Time, &
                          'random overlap cloud amount', &
                          'percent', missing_value=missing_value)

      id_flxnet = &
         register_diag_field (mod_name, 'flxnet', bxes(1:3), Time, &
                          'net longwave radiative flux', &
                          'W/m**2', missing_value=missing_value)

      id_fsw    = &
         register_diag_field (mod_name, 'fsw', bxes(1:3), Time, &
                          'net shortwave radiative flux', &
                          'W/m**2', missing_value=missing_value)

      id_ufsw   = &
         register_diag_field (mod_name, 'ufsw', bxes(1:3), Time, &
                          'upward shortwave radiative flux ', &
                          'W/m**2', missing_value=missing_value)

      id_dfsw   = &
         register_diag_field (mod_name, 'dfsw', bxes(1:3), Time, &
                          'downward shortwave radiative flux ', &
                          'W/m**2', missing_value=missing_value)

      id_psj    = &
         register_diag_field (mod_name, 'psj', axes(1:2), Time, &
                          'surface pressure', &
                          'Pa', missing_value=missing_value)

      id_tmpsfc = &
         register_diag_field (mod_name, 'tmpsfc', axes(1:2), Time, &
                          'surface temperature', &
                          'deg_K', missing_value=missing_value)

      id_cvisrfgd_dir = &
         register_diag_field (mod_name, 'cvisrfgd_dir', axes(1:2), Time , &
                         'direct visible surface albedo', &
                        'dimensionless', missing_value=missing_value)

       id_cvisrfgd_dif = &
       register_diag_field (mod_name, 'cvisrfgd_dif', axes(1:2), Time, &
                          'diffuse visible surface albedo', &
                          'dimensionless', missing_value=missing_value)
 
       id_cirrfgd_dir = &
        register_diag_field (mod_name, 'cirrfgd_dir', axes(1:2), Time, &
                       'direct infra-red surface albedo', &
                       'dimensionless', missing_value=missing_value)

       id_cirrfgd_dif = &
         register_diag_field (mod_name, 'cirrfgd_dif', axes(1:2), Time, &
                       'diffuse infra-red surface albedo', &
                      'dimensionless', missing_value=missing_value)

!-----------------------------------------------------------------------
! Change by Xianglei Huang
!-----------------------------------------------------------------------
       id_lw_bdyflx(1) = &
         register_diag_field (mod_name, 'olr_0_560', axes(1:2), Time, &
                       'olr in 0_560  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx(2) = &
         register_diag_field (mod_name, 'olr_560_800', axes(1:2), Time, &
                       'olr in 560_800  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx(3) = &
         register_diag_field (mod_name, 'olr_800_900', axes(1:2),&
                             Time, 'olr in 800_900 band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx(4) = &
         register_diag_field (mod_name, 'olr_990_1070', axes(1:2), &
                              Time, 'olr in 990_1070 band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx(5) = &
         register_diag_field  (mod_name, 'olr_900_990', axes(1:2), &
                               Time, 'olr in 900_990 band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx(6) = &
         register_diag_field  (mod_name, 'olr_1070_1200', axes(1:2), &
                               Time, 'olr in 1070_1200 band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx(7) = &
         register_diag_field  (mod_name, 'olr_1200_1400', axes(1:2), &
                               Time, 'olr in 1200_1400 band', &
                      'W/m**2', missing_value=missing_value)
      id_bT_Aqua31 = register_diag_field ( mod_name, 'bT_Aqua31', axes(1:2), &
           Time, 'brightness temperature in Aqua 31 885_925 band', 'K', &
           missing_value=missing_value)

!-----------------------------------------------------------------------
! End of Change 
!-----------------------------------------------------------------------

       id_sw_bdyflx(1) = &
         register_diag_field (mod_name, 'swup_toa_vis', axes(1:2),  &
                       Time, &
                       'sw up flx in vis band at toa', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx(2) = &
         register_diag_field (mod_name, 'swup_toa_1p6', axes(1:2),  &
                       Time, &
                       'sw up flx in 1.6 micron band at toa', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx(3) = &
         register_diag_field (mod_name, 'swnt_sfc_vis', axes(1:2),  &
                       Time, &
                       'net sw flx in vis band at sfc', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx(4) = &
         register_diag_field (mod_name, 'swnt_sfc_1p6', axes(1:2),  &
                       Time, &
                       'net sw flx in 1.6 micron band at sfc', &
                      'W/m**2', missing_value=missing_value)

!-------------------------------------------------------------------

      if (do_totcld_forcing) then
        id_radswpcf = &
           register_diag_field (mod_name, 'radswpcf', axes(1:3), Time, &
                            'temperature forcing from sw w/o clouds', &
                            'deg_K/sec', missing_value=missing_value)

        id_radpcf = &
           register_diag_field (mod_name, 'radpcf', axes(1:3), Time, &
                            'temperature forcing w/o clouds', &
                            'deg_K/sec', missing_value=missing_value)

        id_flxnetcf = &
           register_diag_field (mod_name, 'flxnetcf', bxes(1:3), Time, &
                            'net longwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

        id_fswcf = &
           register_diag_field (mod_name, 'fswcf', bxes(1:3), Time, &
                            'net shortwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

        id_ufswcf   = &
           register_diag_field (mod_name, 'ufswcf', bxes(1:3), Time, &
                            'upward shortwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

        id_dfswcf   = &
           register_diag_field (mod_name, 'dfswcf', bxes(1:3), Time, &
                            'downward shortwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

!--------------------------------------------------------------------------
! Change by Xianglei Huang
!--------------------------------------------------------------------------

       id_lw_bdyflx_clr(1) = &
         register_diag_field (mod_name, 'olr_0_560_cf', axes(1:2), Time, &
                       'clr sky olr in 0_560  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(2) = &
         register_diag_field (mod_name, 'olr_560_800_cf', axes(1:2), Time, &
                       'clr sky olr in 560_800  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(3) = &
         register_diag_field (mod_name, 'olr_800_900_cf', axes(1:2),&
                             Time, 'clr sky olr in 800_900 band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(4) = &
         register_diag_field (mod_name, 'olr_990_1070_cf', axes(1:2), &
                              Time, 'clr sky olr in 990_1070 band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(5) = &
         register_diag_field (mod_name, 'olr_900_990_cf', axes(1:2), Time, &
                       'clr sky olr in 900_990  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(6) = &
         register_diag_field (mod_name, 'olr_1070_1200_cf', axes(1:2),&
                             Time, 'clr sky olr in 1070_1200 band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(7) = &
         register_diag_field (mod_name, 'olr_1200_1400_cf', axes(1:2), &
                              Time, 'clr sky olr in 1200_1400 band', &
                      'W/m**2', missing_value=missing_value)

!---------------------------------------------------------------------------
! End of Change 
!---------------------------------------------------------------------------

       id_sw_bdyflx_clr(1) = &
         register_diag_field (mod_name, 'swup_toa_vis_cf', axes(1:2),  &
                       Time, &
                       'clr sky sw up flx in vis band at toa', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx_clr(2) = &
         register_diag_field (mod_name, 'swup_toa_1p6_cf', axes(1:2),  &
                       Time, &
                       'clr sky sw up flx in 1.6 micron band at toa', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx_clr(3) = &
         register_diag_field (mod_name, 'swnt_sfc_vis_cf', axes(1:2),  &
                       Time, &
                       'clr sky net sw flx in vis band at sfc', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx_clr(4) = &
         register_diag_field (mod_name, 'swnt_sfc_1p6_cf', axes(1:2),  &
                       Time, &
                       'clr sky net sw flx in 1.6 micron band at sfc', &
                      'W/m**2', missing_value=missing_value)

      endif

!---------------------------------------------------------------------


end subroutine register_fields


!####################################################################


                  end module rad_output_file_mod


