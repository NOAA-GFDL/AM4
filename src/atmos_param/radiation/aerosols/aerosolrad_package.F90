                 module aerosolrad_package_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    aerosolrad_package_mod provides the radiative properties 
!    associated with the atmospheric aerosols.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
!    shared modules:

use mpp_mod,               only: input_nml_file
use fms_mod,               only: open_namelist_file, fms_init, &
                                 mpp_pe, mpp_root_pe, stdlog, &
                                 file_exist, write_version_number, &
                                 check_nml_error, error_mesg, &
                                 FATAL, NOTE, close_file
use constants_mod,         only: diffac, constants_init
use mpp_io_mod,            only: mpp_open, mpp_close, MPP_RDONLY,   &
                                 MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI, &
                                 MPP_SINGLE, mpp_io_init
use time_manager_mod,      only: time_type, time_manager_init,  &
                                 get_date, set_date, operator(+), &
                                 print_date, operator(-), operator(>)
use diag_manager_mod,      only: diag_manager_init, get_base_time
use interpolator_mod,      only: interpolate_type, interpolator_init, &
                                 interpolator, interpolator_end, &
                                 obtain_interpolator_time_slices, &
                                 unset_interpolator_time_flag, &
                                 CONSTANT, INTERP_WEIGHTED_P
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index, NO_TRACER

! shared radiation package modules:
                                
use aerosol_types_mod,     only: aerosol_type

use aerosolrad_types_mod,  only: aerosolrad_control_type, &
                                 aerosolrad_diag_type

use sealw99_mod,           only: NBLW

use esfsw_driver_mod,      only: esfsw_number_of_bands, &
                                 esfsw_band_segments, &
                                 esfsw_bands, esfsw_thickavg

!-------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    aerosolrad_package_mod provides the radiative properties 
!    associated with the atmospheric aerosols.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public           &
       aerosolrad_package_init, aerosol_radiative_properties, &
       aerosolrad_package_time_vary, aerosolrad_package_endts, &
       aerosolrad_package_end, number_of_lw_aerosol_bands
     
private          &

!  called from aerosolrad_package_init: 
   assign_aerosol_opt_props, read_optical_input_file, &
   sw_aerosol_interaction,  lw_aerosol_interaction
       

!---------------------------------------------------------------------
!-------- namelist  ---------

integer, parameter   ::        &
             MAX_OPTICAL_FIELDS = 1100 ! maximum number of aerosol 
                                       ! optical property types

integer, parameter   ::        &
             NUM_AERO_INDICES   = 12   ! Last dimension of opt_indices
logical              ::        &
             do_lwaerosol = .false.    ! aerosol efects included in lw
                                       ! radiation ?
logical              ::        &
             do_swaerosol = .false.    ! aerosol effects included in sw
                                       ! radiation ?
logical              ::        &    
             force_to_repro_quebec = .false.
                                       ! if true, code sequence is 
                                       ! executed which reproduces
                                       ! quebec+ answers for 
                                       ! AM3p8e in pre_Riga
character(len=48)    ::        &
             aerosol_data_set = ' '    ! source of aerosol data; if 
                                       ! aerosols not desired remains
                                       ! ' ', otherwise is set to either
                                       ! 'shettle_fenn' or 
                                       ! 'Ginoux_Reddy'

!----------------------------------------------------------------------
!    the avaialable aerosol datasets are :
!    1) "shettle_fenn":  
!        Ref: shettle, e.p. and r.w. fenn, models for the aerosols of 
!            the lower atmosphere and the effects of humidity variations
!            on their optical properties,afgl-tr-79-0214,1979,94pp.    
!    2) "Ginoux_Reddy": 3D Aerosol fields are generated online reflecting
!        emissions, transport and deposition:
!        Ref: Reddy et al., 2005, Ginoux et al., 2005
!        
!----------------------------------------------------------------------

character(len=64)    ::        &
             aerosol_optical_names(MAX_OPTICAL_FIELDS) = '  '
                                       ! names associated with the 
                                       ! optical property types that
                                       ! are to be used in this 
                                       ! experiment
character(len=64)    ::        &
             optical_filename = ' '    ! name of file containing the
                                       ! aerosol optical property types
logical              ::        &
             using_volcanic_sw_files = .false.
                                       ! files containing sw aerosol
                                       ! optical properties from vol-
                                       ! canic activity are to be
                                       ! used to supplement those cal-
                                       ! culated by model ?
logical              ::        &
             using_volcanic_lw_files = .false.
                                       ! files containing lw aerosol
                                       ! optical properties from vol-
                                       ! canic activity are to be
                                       ! used to supplement those cal-
                                       ! culated by model ?
character(len=64)    ::        &
              sw_ext_filename = ' '    ! name of file containing the
                                       ! aerosol sw extinction optical
                                       ! depth
character(len=64)    ::        &
              sw_ssa_filename = ' '    ! name of file containing the
                                       ! aerosol sw single scattering 
                                       ! albedo
character(len=64)    ::        &
              sw_asy_filename = ' '    ! name of file containing the
                                       ! aerosol sw asymmetry factor   
character(len=64)    ::        &
              lw_ext_filename = ' '    ! name of file containing the
                                       ! aerosol lw extinction optical
                                       ! depth
character(len=64)    ::        &
              lw_ssa_filename = ' '    ! name of file containing the
                                       ! aerosol lw single scattering 
                                       ! albedo
character(len=64)    ::        &
              lw_asy_filename = ' '    ! name of file containing the
                                       ! aerosol lw asymmetry factor   
                                       ! the supplemental input files
character(len=64)    ::        &
              sw_ext_root                  = '   ' 
                                       ! names given to sw extopdep in
                                       ! input netcdf file
character(len=64)    ::        &
              sw_ssa_root                  = '   ' 
                                       ! name given to sw single scat-
                                       ! tering albedo in input netcdf 
                                       ! file
character(len=64)    ::        &
              sw_asy_root                  = '   ' 
                                       ! name given to sw asymmetry
                                       ! factor in input netcdf file
character(len=64)    ::        &
              lw_ext_root                  = '   ' 
                                       ! name given to lw extopdep in
                                       ! input netcdf file
character(len=64)    ::        &
              lw_ssa_root                  = '   '  
                                       ! name given to lw single scat-
                                       ! tering albedo in input netcdf 
                                       ! file
character(len=64)    ::        &
              lw_asy_root                  = '   '      
                                       ! name given to lw asymmetry
                                       ! factor in input netcdf file
integer, dimension(6) ::       &
              volcanic_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /) 
                                       ! time in volcanic data set
                                       ! corresponding to model
                                       ! initial time 
                                       ! (yr, mo, dy, hr, mn, sc)
logical :: interpolating_volcanic_data = .true.
                                       ! volcanic datasets will be
                                       ! time interpolated rather than
                                       ! held constant for a month ?
logical :: repeat_volcano_year = .false. 
                                      ! the same single year's data from
                                      ! the input data set should be 
                                      ! used for each model year ?
integer :: volcano_year_used = 0      ! year of volcanic data to repeat
                                      ! when repeat_volcano_year is
                                      ! .true.
logical :: using_im_bcsul = .false.   ! bc and sulfate aerosols are 
                                      ! treated as an internal mixture ?
integer, dimension(0:100) ::  nitrate_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  omphilic_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  bcphilic_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt1_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt2_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt3_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt4_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt5_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt_aitken_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt_fine_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt_coarse_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  sulfate_indices = (/      &
                           30,30,30,30,30,30,30,30,30,30,30,30,30,30, &
                           30,30,30,30,30,30,30,30,30,30,30,30,30,30, &
                           30,30,30,30,30,35,35,35,35,35,40,40,40,40, &
                           40,45,45,45,45,45,50,50,50,50,50,55,55,55, &
                           55,55,60,60,60,60,60,65,65,65,65,65,70,70, &
                           70,70,70,75,75,75,75,75,80,80,80,80,82,82, &
                           84,84,86,86,88,88,90,91,92,93,94,95,96,97, &
                           98,99,100 /)
!yim
integer, dimension(0:100) ::  sulfate_vol_indices = (/      &
                             100,98,98,96,96,94,94,92,92,90,90,88,88,86,86,84,84,82,82,80, &
                             80,80,80,75,75,75,75,75,70,70,70,70,70,65,65,65, &
                             65,65,60,60,60,60,60,55,55,55,55,55,50,50,50,50,50,45,45,45, &
                             45,45,40,40,40,40,40,35,35,35,35,35,30,30,30,30,30,25,25,25, &
                             25,25,20,20,20,20,20,15,15,15,15,15,10,10,10,10,10,5,5,5, &
                             5,5,0,0,0  /)


namelist / aerosolrad_package_nml /                          &
                                    do_lwaerosol, do_swaerosol, &
                                    force_to_repro_quebec, &
                                    aerosol_data_set, &
                                    aerosol_optical_names, &
                                    sulfate_indices, &
                                    sulfate_vol_indices, &
                                    omphilic_indices, &
                                    bcphilic_indices, &
                                    seasalt1_indices, &
                                    seasalt2_indices, &
                                    seasalt3_indices, &
                                    seasalt4_indices, &
                                    seasalt5_indices, &
                                    seasalt_aitken_indices, &
                                    seasalt_fine_indices, &
                                    seasalt_coarse_indices, &
                                    nitrate_indices,  &
                                    optical_filename   , &
                                    using_volcanic_sw_files, &
                                    using_volcanic_lw_files, &
                                    volcanic_dataset_entry, &
                                    interpolating_volcanic_data, &
                                    repeat_volcano_year, &
                                    volcano_year_used, &
                                    using_im_bcsul, &
                                    sw_ext_filename, sw_ssa_filename, &
                                    sw_asy_filename, lw_ext_filename, &
                                    lw_ssa_filename, lw_asy_filename, &
                                    sw_ext_root, sw_ssa_root,   &
                                    sw_asy_root, lw_ext_root,   &
                                    lw_ssa_root, lw_asy_root

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


!---------------------------------------------------------------------
!   the following are interpolate_type variables containing the
!   additional aerosol optical properties that may be included as
!   input to the radiation package.
!---------------------------------------------------------------------
type(interpolate_type), save  :: Sw_aer_extopdep_interp
type(interpolate_type), save  :: Sw_aer_ssalb_interp
type(interpolate_type), save  :: Sw_aer_asymm_interp
type(interpolate_type), save  :: Lw_aer_extopdep_interp
type(interpolate_type), save  :: Lw_aer_ssalb_interp
type(interpolate_type), save  :: Lw_aer_asymm_interp

!---------------------------------------------------------------------
!    the following variables define the number and type of different
!    bands over which the radiation package calculates aerosol 
!    radiative properties.
!---------------------------------------------------------------------
integer, parameter ::    &
         N_AEROSOL_BANDS_FR = 8 ! number of non-continuum ir aerosol
                                ! emissivity bands 
integer, parameter ::     &
         N_AEROSOL_BANDS_CO = 1 ! number of continuum ir aerosol
                                ! emissivity bands  
integer, parameter ::     &
         N_AEROSOL_BANDS_CN = 1 ! number of diagnostic continuum ir 
                                ! aerosol emissivity bands  
integer, parameter ::    &
         N_AEROSOL_BANDS = N_AEROSOL_BANDS_FR + N_AEROSOL_BANDS_CO
                                ! total number of ir aerosol emissivity
                                ! bands 

!--------------------------------------------------------------------
!    num_wavenumbers is the number of wavenumber bands over which the
!    aerosol parameterization provides aerosol radiative property data.
!--------------------------------------------------------------------
integer     ::  num_wavenumbers = 0 ! number of wavenumber bands 
                                    ! present in the aerosol 
                                    ! parameterization

!----------------------------------------------------------------------
!    the following variable defines the number of aerosol property 
!    types that are active.
!----------------------------------------------------------------------
integer     :: naermodels = 0   ! number of aerosol optical properties
                                ! types that are active

!---------------------------------------------------------------------
!    flags indicating an index value characteristic of the optical prop-
!    erties associated with different aerosols
!---------------------------------------------------------------------
integer, PARAMETER ::   SULFATE_FLAG =  0
integer, PARAMETER ::  OMPHILIC_FLAG = -1
integer, PARAMETER ::  BCPHILIC_FLAG = -2
integer, PARAMETER ::  SEASALT1_FLAG = -3
integer, PARAMETER ::  SEASALT2_FLAG = -4
integer, PARAMETER ::  SEASALT3_FLAG = -5
integer, PARAMETER ::  SEASALT4_FLAG = -6
integer, PARAMETER ::  SEASALT5_FLAG = -7
integer, PARAMETER ::  SEASALTA_FLAG = -8
integer, PARAMETER ::  SEASALTF_FLAG = -9
integer, PARAMETER ::  SEASALTC_FLAG = -10
integer, PARAMETER ::  NITRATE_FLAG  = -11
! Add additional aerosols before BC_FLAG which is a dummy. 
! Then BC_FLAG should be the lowest flag number-1
!yim
integer, PARAMETER ::  BC_FLAG = -12 
!integer, PARAMETER ::  NOT_IN_USE = -2000

!----------------------------------------------------------------------
!    the following index arrays contain the mapping information between 
!    actual model relative humidity and the available enties in the 
!    aerosol optical properties file.
!----------------------------------------------------------------------
integer, dimension(:,:), allocatable :: sulfate_index_MOD
integer, dimension(:),   allocatable :: optical_index_MOD
integer, dimension(:),   allocatable :: omphilic_index_MOD
integer, dimension(:),   allocatable :: bcphilic_index_MOD
integer, dimension(:),   allocatable :: seasalt1_index_MOD
integer, dimension(:),   allocatable :: seasalt2_index_MOD
integer, dimension(:),   allocatable :: seasalt3_index_MOD
integer, dimension(:),   allocatable :: seasalt4_index_MOD
integer, dimension(:),   allocatable :: seasalt5_index_MOD
integer, dimension(:),   allocatable :: seasalt_aitken_index_MOD
integer, dimension(:),   allocatable :: seasalt_fine_index_MOD
integer, dimension(:),   allocatable :: seasalt_coarse_index_MOD
integer, dimension(:),   allocatable :: nitrate_index_MOD

!---------------------------------------------------------------------
!    the following arrays related to sw aerosol effects are allocated 
!    during initialization and retained throughout the integration.
!    here n refers to the bands of the solar parameterization, ni
!    to the bands of the aerosol parameterization, and na to the optical
!    properties type.
!
!      solivlaero(n,ni)  amount of toa incoming solar from solar
!                        spectral band n that is in aerosol parameter-
!                        ization band ni
!      nivl1aero(n)      the aerosol band index corresponding to the 
!                        lowest wave number of spectral band n
!      nivl2aero(n)      the aerosol band index corresponding to the 
!                        highest wave number of spectral band n
!      endaerwvnsf(ni)   ending wave number of aerosol parameterization
!                        band ni
!      aeroextivl(ni,na) extinction coefficient for aerosol parameter-
!                        ization band ni for aerosol optical property 
!                        type na
!      aerossalbivl(ni,na) 
!                        single-scattering albedo for aerosol band 
!                        ni and aerosol optical property type na
!      aeroasymmivl(ni,na)
!                        asymmetry factor for aerosol band ni  and 
!                        aerosol optical property type na
!
!---------------------------------------------------------------------
real,    dimension(:,:), allocatable   :: solivlaero  
integer, dimension(:),   allocatable   :: nivl1aero, nivl2aero
integer, dimension(:),   allocatable   :: endaerwvnsf
real,    dimension(:,:), allocatable   :: aeroextivl, aerossalbivl, &
                                          aeroasymmivl

!---------------------------------------------------------------------
!    sfl following arrays related to lw aerosol effects are allocated 
!    during initialization and retained throughout the integration.
!
!    sflwwts(n,ni)     the fraction of the planck function in aerosol 
!                      emissivity band n that is in aerosol param-
!                      eterization band ni
!
!----------------------------------------------------------------------
real,    dimension(:,:), allocatable   :: sflwwts, sflwwts_cn

!--------------------------------------------------------------------
!    logical flags 
!--------------------------------------------------------------------
logical :: module_is_initialized      = .false. ! module has been
                                                ! initialized ?
!logical :: doing_predicted_aerosols   = .false. ! predicted aerosol 
                                                ! scheme being used ?
logical :: band_calculation_completed = .false. ! lw properties have
                                                ! been calculated ?

type(time_type) :: Volcanic_offset  ! difference between model initial
                                    ! time and volcanic timeseries app-
                                    ! lied at model initial time
                                    ! [ time_type ]
type(time_type) :: Volcanic_entry   ! time in volcanic timeseries which
                                    ! is mapped to model initial time
                                    ! [ time_type ]
logical    :: negative_offset = .false.
                                !  the model initial time is later than
                                !  the volcanic_dataset_entry time  ?
integer :: nfields_sw_ext = 0   ! number of fields contained in 
                                ! supplemental sw_ext file
integer :: nfields_sw_ssa = 0   ! number of fields contained in 
                                ! supplemental sw_ssa file
integer :: nfields_sw_asy = 0   ! number of fields contained in 
                                ! supplemental sw_asy file
integer :: nfields_lw_ext = 0   ! number of fields contained in 
                                ! supplemental lw_ext file
integer :: nfields_lw_ssa = 0   ! number of fields contained in 
                                ! supplemental lw_ssa file
integer :: nfields_lw_asy = 0   ! number of fields contained in 
                                ! supplemental lw_asy file

!-------------------------------------------------------------------
!   arrays holding variable names:
character(len=64), dimension(:), allocatable ::   &
                                sw_ext_name, sw_ssa_name, sw_asy_name, &
                                lw_ext_name, lw_ssa_name, lw_asy_name

!-------------------------------------------------------------------
!    arrays to hold data when not interpolating on every step:
real, dimension(:,:,:,:), allocatable :: sw_ext_save
real, dimension(:,:,:,:), allocatable :: sw_ssa_save
real, dimension(:,:,:,:), allocatable :: sw_asy_save
real, dimension(:,:,:,:), allocatable :: lw_ext_save
real, dimension(:,:,:,:), allocatable :: lw_ssa_save
real, dimension(:,:,:,:), allocatable :: lw_asy_save

!---------------------------------------------------------------------
!   module variables to hold values unchanging in time:
!---------------------------------------------------------------------
real, dimension(:,:), allocatable :: aerextband_MOD
real, dimension(:,:), allocatable :: aerssalbband_MOD
real, dimension(:,:), allocatable :: aerasymmband_MOD
real, dimension(:,:), allocatable :: aerextbandlw_MOD
real, dimension(:,:), allocatable :: aerssalbbandlw_MOD
real, dimension(:,:), allocatable :: aerextbandlw_cn_MOD
real, dimension(:,:), allocatable :: aerssalbbandlw_cn_MOD

!---------------------------------------------------------------------
!    logical variables indicating whether interpolation is currently
!    needed:
logical :: need_sw_ext = .true.
logical :: need_sw_ssa = .true.
logical :: need_sw_asy = .true.
logical :: need_lw_ext = .true.
logical :: need_lw_ssa = .true.
logical :: need_lw_asy = .true.

!---------------------------------------------------------------------
!    logical variables indicating whether the particular radiative 
!    property associated with volcanoes is being supplied:
logical :: using_sw_ext = .false. 
logical :: using_sw_ssa = .false. 
logical :: using_sw_asy = .false. 
logical :: using_lw_ext = .false. 
logical :: using_lw_ssa = .false. 
logical :: using_lw_asy = .false. 

!---------------------------------------------------------------------
!    counters associated with determining when interpolation needs to
!    be done:
logical :: mo_save_set = .false.
integer :: mo_save = 0
integer :: mo_new

type(time_type) :: Volcano_time
integer :: nfields_save
integer :: num_sul, num_bc
integer, dimension(:), allocatable :: sul_ind, bc_ind

!---------------------------------------------------------------------
 


                         contains
 

 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="aerosolrad_package_init">
!  <OVERVIEW>
!     aerosolrad_package_init is the constructor for 
!     aerosolrad_package_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!     aerosolrad_package_init is the constructor for 
!     aerosolrad_package_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosolrad_package_init (aerosol_names)
!  </TEMPLATE>
!  <IN NAME="aerosol_names" TYPE="character">
!   names of the activated aerosol species
!  </IN>
! </SUBROUTINE>
!
subroutine aerosolrad_package_init (kmax, aerosol_names, lonb, latb, &
                                    Aerosolrad_control)

!---------------------------------------------------------------------
!     aerosolrad_package_init is the constructor for 
!     aerosolrad_package_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!character(len=64), dimension(:), intent(in)  :: aerosol_names
integer,                        intent(in)    :: kmax
character(len=*), dimension(:), intent(in)    :: aerosol_names
real, dimension(:,:),           intent(in)    :: lonb,latb
type(aerosolrad_control_type),  intent(inout) :: Aerosolrad_control

!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      kmax              number of model levels
!      aerosol_names     the names assigned to each of the activated
!                        aerosol species
!       lonb           2d array of model longitudes at cell corners
!                      [ radians ]
!       latb           2d array of model latitudes at cell corners
!                      [ radians ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer        :: unit, ierr, io, logunit
      integer        :: n
      character(len=16) :: chvers
      character(len=4)  :: chyr  
      integer           :: nbands    ! number of solar bands

      type(time_type) :: Time_init  ! initial calendar time for model  
                                    ! [ time_type ]

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
      call mpp_io_init
      call fms_init
      call constants_init
      call diag_manager_init
      call time_manager_init

       nfields_save = size(aerosol_names(:))

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=aerosolrad_package_nml, iostat=io)
      ierr = check_nml_error(io,'aerosolrad_package_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosolrad_package_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'aerosolrad_package_nml')
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
                          write (logunit, nml=aerosolrad_package_nml)

!----------------------------------------------------------------------
!    define control variables which indicate whether the impact of 
!    aerosols on radiation is to be included in the sw and lw rad-
!    iation calculations. define a control variable which will be true 
!    if aerosols are included in either the sw or the lw radiation 
!    (Aerosolrad_control%do_aerosol).
!----------------------------------------------------------------------
      Aerosolrad_control%do_swaerosol = do_swaerosol
      Aerosolrad_control%do_lwaerosol = do_lwaerosol
      if (do_lwaerosol .or. do_swaerosol  .or.  &
          using_volcanic_lw_files .or. using_volcanic_sw_files .or. &
         Aerosolrad_control%do_lwaerosol_forcing .or.  &
         Aerosolrad_control%do_swaerosol_forcing)  then
        Aerosolrad_control%do_aerosol = .true.
      else
        Aerosolrad_control%do_aerosol = .false.
      endif

!--------------------------------------------------------------------
!     number of solar bands
!---------------------------------------------------------------------
      call esfsw_number_of_bands (nbands)
!---------------------------------------------------------------------
!    exit if an aerosol_data_set is provided when do_aerosol is 
!    .false..
!---------------------------------------------------------------------
      if ( .not. Aerosolrad_control%do_aerosol .and.    &
           trim(aerosol_data_set) /= ' ') then
        call error_mesg ('aerosolrad_package_mod', &
           'if aerosol impacts are not desired, aerosol_data_set '//&
            'must be set to "   "', FATAL)
      endif

!---------------------------------------------------------------------
!    exit if no aerosol_data_set is provided when do_aerosol  
!    is .true..
!---------------------------------------------------------------------
      if ( Aerosolrad_control%do_aerosol .and.    &
           trim(aerosol_data_set) == ' ') then
        call error_mesg ('aerosolrad_package_mod', &
           'if aerosol impacts are desired, aerosol_data_set '//&
            'must be non-blank', FATAL)
      endif

!---------------------------------------------------------------------
!    exit if aerosol effects are desired but the aerosol input file
!    provided no aerosol fields.
!---------------------------------------------------------------------
      if (Aerosolrad_control%do_aerosol .and. size(aerosol_names(:)) == 0) then
        call error_mesg ('aerosolrad_package_mod', &
          ' aerosols desired  for radiation but no aerosol '//&
            'data_names supplied', FATAL)
      endif


!----------------------------------------------------------------------
!    if aerosol radiative effects are to be included, call 
!    assign_aerosol_opt_props to assign the proper aerosol 
!    properties type to each aerosol type. then call 
!    read_optical_input_file to read the optical input file contain-
!    ing the aerosol parameterization information and data.
!----------------------------------------------------------------------
      if (Aerosolrad_control%do_aerosol) then
        call assign_aerosol_opt_props (aerosol_names)
        call read_optical_input_file
      endif
 
!---------------------------------------------------------------------
!    if aerosol effects are to be included in the sw calculation,
!    call sw_aerosol_interaction to define the weights needed to prop-
!    erly map the input data from the aerosol parameterization bands to 
!    the solar parameterization bands that the model is using.
!--------------------------------------------------------------------
      if (do_swaerosol .or. Aerosolrad_control%do_swaerosol_forcing) then
        call sw_aerosol_interaction                 
      endif

!---------------------------------------------------------------------
!    if aerosol effects are to be included in the lw calculation,
!    call lw_aerosol_interaction to define the weights needed to prop-
!    erly map the input data from the aerosol parameterization bands to 
!    the solar parameterization bands that the model is using. if
!    they are not, indicate that this part of the code has been 
!    executed.
!---------------------------------------------------------------------
      if (do_lwaerosol .or. Aerosolrad_control%do_lwaerosol_forcing) then
        call lw_aerosol_interaction
      endif

!---------------------------------------------------------------------
!    make sure consistent nml settings are present. Cannot use volcanic
!    aerosols unless model aerosols are also activated.
!---------------------------------------------------------------------
      if (.not. do_swaerosol  .and.   &
          using_volcanic_sw_files) then
        call error_mesg ('aerosolrad_package_mod', &
         'cant use sw volcanic aerosols without activating standard &
                                               & sw aerosols', FATAL)
      endif
      if (.not. do_lwaerosol  .and.   &
          using_volcanic_lw_files) then
        call error_mesg ('aerosolrad_package_mod', &
         'cant use lw volcanic aerosols without activating standard &
                                               & lw aerosols', FATAL)
      endif

!---------------------------------------------------------------------
!    set the volcanic control variables to .false. when the model 
!    aerosols are not active.
!---------------------------------------------------------------------
      Aerosolrad_control%volcanic_lw_aerosols = using_volcanic_lw_files
      Aerosolrad_control%volcanic_sw_aerosols = using_volcanic_sw_files

!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the volcanic files are
!    to be used; if not present, the dataset entry point is taken as
!    the model base_time, defined in the diag_table.
!---------------------------------------------------------------------
      Time_init = get_base_time()
      if (using_volcanic_sw_files .or.  &
          using_volcanic_lw_files) then
        if (volcanic_dataset_entry(1) == 1 .and. &
            volcanic_dataset_entry(2) == 1 .and. &
            volcanic_dataset_entry(3) == 1 .and. &
            volcanic_dataset_entry(4) == 0 .and. &
            volcanic_dataset_entry(5) == 0 .and. &
            volcanic_dataset_entry(6) == 0 ) then      
          Volcanic_entry = Time_init

!----------------------------------------------------------------------
!    define the offset from model base time  (defined in diag_table) 
!    to volcanic_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
      else
        Volcanic_entry  = set_date (volcanic_dataset_entry(1), &
                                    volcanic_dataset_entry(2), &
                                    volcanic_dataset_entry(3), &
                                    volcanic_dataset_entry(4), &
                                    volcanic_dataset_entry(5), &
                                    volcanic_dataset_entry(6))
       endif
     else
       Volcanic_entry  = set_date (volcanic_dataset_entry(1), &
                                   volcanic_dataset_entry(2), &
                                   volcanic_dataset_entry(3), &
                                   volcanic_dataset_entry(4), &
                                   volcanic_dataset_entry(5), &
                                   volcanic_dataset_entry(6))

      endif
      if (using_volcanic_sw_files .or.  &
          using_volcanic_lw_files) then
        if (repeat_volcano_year) then
          if (volcano_year_used == 0) then
            call error_mesg ('aerosolrad_package_init', &
              'valid year must be supplied when &
                       &repeat_volcano_year is .true.', FATAL)
          endif
        endif
        call print_date(Volcanic_entry , str='Data from volcano &
                                           &timeseries at time:')
        call print_date(Time_init , str='This data is mapped to &
                                                  &model time:')
        if (repeat_volcano_year) then
          write (chyr, '(i4)') volcano_year_used
          call error_mesg ('aerosolrad_package_init', &
           'volcanic data from dataset year '  // chyr // ' will be &
                                    &used for all model years.', NOTE)
        endif
      endif
      Volcanic_offset = Volcanic_entry - Time_init

      if (Time_init > Volcanic_entry) then
        negative_offset = .true.
      else
        negative_offset = .false.
      endif

!-----------------------------------------------------------------------
!    if desired, process the sw extinction coefficient file. allocate 
!    space for and define the names of each variable. if not interpol-
!    ating the data, allocate an array to store it between timesteps. 
!    call interpolator_init to initialize the interpolation module for
!    the file.
!-----------------------------------------------------------------------
      if (using_volcanic_sw_files) then
        if (trim(sw_ext_root) /= ' '  ) then
          nfields_sw_ext = nbands ! num of solar bands
          allocate (sw_ext_name (nfields_sw_ext))
          if (.not. interpolating_volcanic_data) then
            allocate (sw_ext_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_sw_ext) )
            sw_ext_save = 0.0
          endif
          do n=1, nfields_sw_ext
            if (n<= 9) then
              write (chvers, '(i1)') n
             sw_ext_name(n) = trim(sw_ext_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              sw_ext_name(n) = trim(sw_ext_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Sw_aer_extopdep_interp,  &
                                  sw_ext_filename, lonb, latb,  &
                                  sw_ext_name(:nfields_sw_ext),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) )  
          using_sw_ext = .true.
        endif

!--------------------------------------------------------------------
!    if desired, process the sw single scattering albedo file. allocate 
!    space for and define the names of each variable. if not interpol-
!    ating the data, allocate an array to store it between timesteps. 
!    call interpolator_init to initialize the interpolation module for
!    the file.
!-----------------------------------------------------------------------
        if (trim(sw_ssa_root) /= ' '  ) then
          nfields_sw_ssa = nbands ! num of solar bands
          allocate (sw_ssa_name (nfields_sw_ssa))
          if (.not. interpolating_volcanic_data) then
            allocate (sw_ssa_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_sw_ssa) )
            sw_ssa_save = 0.0
          endif
          do n=1, nfields_sw_ssa
            if (n<= 9) then
              write (chvers, '(i1)') n
             sw_ssa_name(n) = trim(sw_ssa_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              sw_ssa_name(n) = trim(sw_ssa_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Sw_aer_ssalb_interp,   &
                                  sw_ssa_filename,  lonb, latb,    &
                                  sw_ssa_name(:nfields_sw_ssa),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) ) 
          using_sw_ssa = .true.
        endif

!--------------------------------------------------------------------
!    if desired, process the sw asymmetry factor file. allocate 
!    space for and define the names of each variable. if not interpol-
!    ating the data, allocate an array to store it between timesteps. 
!    call interpolator_init to initialize the interpolation module for
!    the file.
!-----------------------------------------------------------------------
        if (trim(sw_asy_root)    /= ' '  ) then
          nfields_sw_asy = nbands ! num of solar bands
          allocate (sw_asy_name (nfields_sw_asy))
          if (.not. interpolating_volcanic_data) then
            allocate (sw_asy_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_sw_asy) )
            sw_asy_save = 0.0
          endif
          do n=1, nfields_sw_asy
            if (n<= 9) then
              write (chvers, '(i1)') n
             sw_asy_name(n) = trim(sw_asy_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              sw_asy_name(n) = trim(sw_asy_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Sw_aer_asymm_interp,   &
                                  sw_asy_filename, lonb, latb,   &
                                  sw_asy_name(:nfields_sw_asy),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) )  
          using_sw_asy = .true.
        endif
      endif

!-----------------------------------------------------------------------
!    if desired, process the lw extinction coefficient file. allocate 
!    space for and define the names of each variable. if not interpol-
!    ating the data, allocate an array to store it between timesteps. 
!    call interpolator_init to initialize the interpolation module for
!    the file.
!-----------------------------------------------------------------------
      if (using_volcanic_lw_files) then
        if (trim(lw_ext_root)    /= ' '  ) then
          nfields_lw_ext = N_AEROSOL_BANDS
          allocate (lw_ext_name (nfields_lw_ext))
          if (.not. interpolating_volcanic_data) then
            allocate (lw_ext_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_lw_ext) )
            lw_ext_save = 0.0
          endif
          do n=1, nfields_lw_ext
            if (n<= 9) then
              write (chvers, '(i1)') n
             lw_ext_name(n) = trim(lw_ext_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              lw_ext_name(n) = trim(lw_ext_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Lw_aer_extopdep_interp,   &
                                  lw_ext_filename, lonb,  latb,  &
                                  lw_ext_name(:nfields_lw_ext),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) )
          using_lw_ext= .true.
        endif

!--------------------------------------------------------------------
!    if desired, process the lw single scattering albedo file.  it 
!    currently is not needed with the sea lw radiation package. allocate
!    space for and define the names of each variable. call 
!    interpolator_init to initialize the interpolation module for the 
!    file.
!-----------------------------------------------------------------------
        if (trim(lw_ssa_root)    /= ' '  ) then
          nfields_lw_ssa = N_AEROSOL_BANDS
          allocate (lw_ssa_name (nfields_lw_ssa))
          if (.not. interpolating_volcanic_data) then
            allocate (lw_ssa_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_lw_ssa) )
            lw_ssa_save = 0.0
          endif
          do n=1, nfields_lw_ssa
            if (n<= 9) then
              write (chvers, '(i1)') n
             lw_ssa_name(n) = trim(lw_ssa_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              lw_ssa_name(n) = trim(lw_ssa_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Lw_aer_ssalb_interp,  &
                                  lw_ssa_filename, lonb, latb,  &
                                  lw_ssa_name(:nfields_lw_ssa),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) )  
          using_lw_ssa = .true.
        endif

!--------------------------------------------------------------------
!    if desired, process the lw asymmetry factor file.  it currently is
!    not needed with the sea lw radiation package. allocate space for 
!    and define the names of each variable. call interpolator_init to
!    initialize the interpolation module for the file.
!-----------------------------------------------------------------------
        if (trim(lw_asy_root)    /= ' '  ) then
          nfields_lw_asy = N_AEROSOL_BANDS
          allocate (lw_asy_name (nfields_lw_asy))
          if (.not. interpolating_volcanic_data) then
            allocate (lw_asy_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_lw_asy) )
            lw_asy_save = 0.0
          endif
          do n=1, nfields_lw_asy
            if (n<= 9) then
              write (chvers, '(i1)') n
             lw_asy_name(n) = trim(lw_asy_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              lw_asy_name(n) = trim(lw_asy_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Lw_aer_asymm_interp,   &
                                   lw_asy_filename,  lonb, latb,  &
                                   lw_asy_name(:nfields_lw_asy),   &
                                   data_out_of_bounds=(/CONSTANT/), &
                                   vert_interp=(/INTERP_WEIGHTED_P/) )  
          using_lw_asy = .true.
        endif
     endif

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------


end subroutine aerosolrad_package_init


!##########################################################################

subroutine aerosolrad_package_time_vary (Time, Aerosolrad_control)

!-------------------------------------------------------------------------
!    aerosolrad_package_time_vary performs time-dependent, space-independent
!    caluclations for this module
!-------------------------------------------------------------------------

type(time_type),                intent(in) :: Time
type(aerosolrad_control_type),  intent(in) :: Aerosolrad_control


      integer  :: yr, mo, dy, hr, mn, sc
      integer  :: na, ni, nw
    
!---------------------------------------------------------------------
!    define the time for which the volcanic properties will be obtained.
!---------------------------------------------------------------------
        if (using_volcanic_sw_files .or.   &
            using_volcanic_lw_files) then
          if (negative_offset) then
             Volcano_time = Time - Volcanic_offset
          else 
             Volcano_time = Time + Volcanic_offset
          endif
          if (repeat_volcano_year) then
            call get_date (Volcano_time, yr, mo, dy, hr, mn, sc)
            Volcano_time = set_date (volcano_year_used, mo,dy,hr,mn,sc)
          endif

!--------------------------------------------------------------------
!    decide whether the volcanic data must be interpolated on this step.
!    if interpolating_volcanic_data is true, then all variables will
!    always be interpolated. when this is not .true., determine if the
!    month of the data desired has changed from the previous value. if
!    it has set the Volcano_time to 12Z on the 15th of the month, and
!    indicate that new data is needed. On the initial call of the job,
!    one always obtains the data (mo_save_set = .false.).
!--------------------------------------------------------------------
          if (interpolating_volcanic_data) then
            need_sw_ext = .true.
            need_sw_ssa = .true.
            need_sw_asy = .true.
            need_lw_ext = .true.
            need_lw_ssa = .true.
            need_lw_asy = .true.
          else
            call get_date (Volcano_time, yr,mo,dy,hr,mn,sc)
            Volcano_time =  set_date (yr, mo,15,12,0,0)
            if (mo_save_set) then
              if (mo /= mo_save) then
                mo_new = mo       
                need_sw_ext = .true.
                need_sw_ssa = .true.
                need_sw_asy = .true.
                need_lw_ext = .true.
                need_lw_ssa = .true.
                need_lw_asy = .true.
              endif
            else
              need_sw_ext = .true.
              need_sw_ssa = .true.
              need_sw_asy = .true.
              need_lw_ext = .true.
              need_lw_ssa = .true.
              need_lw_asy = .true.
            endif
          endif
        endif ! (using_volcanic_lw or using_volcanic_sw)

!--------------------------------------------------------------------
!    if the volcanic sw aerosol extinction is being supplied, make sure
!    needed time slices are available.
!--------------------------------------------------------------------
        if (using_sw_ext) then
          if (need_sw_ext) then
            if (nfields_sw_ext >= 1) then
              call obtain_interpolator_time_slices  &
                              (Sw_aer_extopdep_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic sw aerosol single scattering albedo is being 
!    supplied, make sure needed time slices are available.
!--------------------------------------------------------------------
        if (using_sw_ssa) then
          if (need_sw_ssa) then
            if (nfields_sw_ssa >= 1) then
              call obtain_interpolator_time_slices  &
                              (Sw_aer_ssalb_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic sw aerosol asymmetry factor is being supplied, 
!    make sure needed time slices are available.
!--------------------------------------------------------------------
        if (using_sw_asy) then
          if (need_sw_asy) then
            if (nfields_sw_asy >= 1) then
              call obtain_interpolator_time_slices  &
                              (Sw_aer_asymm_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw aerosol extinction is being supplied, 
!    make sure needed time slices are available.
!--------------------------------------------------------------------
        if (using_lw_ext) then
          if (need_lw_ext) then
            if (nfields_lw_ext >= 1) then
              call obtain_interpolator_time_slices  &
                              (Lw_aer_extopdep_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw single scattering albedo is being supplied, 
!    make sure needed time slices are available.
!--------------------------------------------------------------------
        if (using_lw_ssa) then
          if (need_lw_ssa) then
            if (nfields_lw_ssa >= 1) then
              call obtain_interpolator_time_slices  &
                              (Lw_aer_ssalb_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw aerosol asymmetry factor is being supplied, 
!    obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_lw_asy) then
          if (need_lw_asy) then
            if (nfields_lw_asy >= 1) then
              call obtain_interpolator_time_slices  &
                              (Lw_aer_asymm_interp, Volcano_Time)
            endif
          endif
        endif

        if (do_lwaerosol .or. Aerosolrad_control%do_lwaerosol_forcing) then
          if (force_to_repro_quebec) then
            if (.not. band_calculation_completed) then
              aerextbandlw_MOD = 0.0
              aerssalbbandlw_MOD = 0.0
              aerextbandlw_cn_MOD = 0.0
              aerssalbbandlw_cn_MOD = 0.0
              do nw=1,naermodels
                do na=1,N_AEROSOL_BANDS
                  do ni=1,num_wavenumbers
                    aerextbandlw_MOD(na,nw) = aerextbandlw_MOD(na,nw) + &
                                aeroextivl(ni,nw)*sflwwts(na,ni)*1.0E+03
                    aerssalbbandlw_MOD(na,nw) = aerssalbbandlw_MOD(na,nw) +   &
                                     aerossalbivl(ni,nw)*sflwwts(na,ni)
                  end do
                end do
              end do
              do nw=1,naermodels
                do na=1,N_AEROSOL_BANDS_CN
                  do ni=1,num_wavenumbers
                    aerextbandlw_cn_MOD(na,nw) = aerextbandlw_cn_MOD(na,nw) + &
                            aeroextivl(ni,nw)*sflwwts_cn(na,ni)*1.0E+03
                    aerssalbbandlw_cn_MOD(na,nw) = aerssalbbandlw_cn_MOD(na,nw) +&
                                  aerossalbivl(ni,nw)*sflwwts_cn(na,ni)
                  end do
                end do
              end do
              band_calculation_completed = .true.
            endif
          endif
        endif

!---------------------------------------------------------------------------

end subroutine aerosolrad_package_time_vary 

!####################################################################

subroutine aerosolrad_package_endts

!---------------------------------------------------------------------
!    when data is not always interpolated, set flags to indicate whether
!    data must be obtained on the next call to this subroutine. if 
!    the current call has obtained data, set the flag indicating that 
!    data is not needed on the next call. also set the flag to indicate 
!    that the initial call has been completed (mo_save_set), and that 
!    the month for which data was obtained has been defined (mo_save). 
!---------------------------------------------------------------------
        if (.not. interpolating_volcanic_data) then
          if (need_sw_ext) then
              need_sw_ext = .false.
              need_sw_ssa = .false.
              need_sw_asy = .false.
              need_lw_ext = .false.
              mo_save_set = .true.
              mo_save = mo_new
          endif
        endif

      if (using_sw_ext) then
        call unset_interpolator_time_flag (Sw_aer_extopdep_interp)
      endif
      if (using_sw_ssa) then
        call unset_interpolator_time_flag (Sw_aer_ssalb_interp)
      endif
      if (using_sw_asy) then
        call unset_interpolator_time_flag (Sw_aer_asymm_interp)
      endif
      if (using_lw_ext) then
        call unset_interpolator_time_flag (Lw_aer_extopdep_interp)
      endif
      if (using_sw_ssa) then
        call unset_interpolator_time_flag (Lw_aer_ssalb_interp)
      endif
      if (using_sw_asy) then
        call unset_interpolator_time_flag (Lw_aer_asymm_interp)
      endif

end subroutine aerosolrad_package_endts


!####################################################################
! <SUBROUTINE NAME="aerosol_radiative_properties">
!  <OVERVIEW>
!    aerosol_radiative_properties defines and returns the radiative
!    properties for each aerosol properties type and for each solar
!    parameterization band in the shortwave and for each aerosol 
!    emissivity band in the longwave.
!  </OVERVIEW>
!  <DESCRIPTION>
!    aerosol_radiative_properties defines and returns the radiative
!    properties for each aerosol properties type and for each solar
!    parameterization band in the shortwave and for each aerosol 
!    emissivity band in the longwave.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_radiative_properties (is, ie, js, je, &
!                                         Aerosolrad_diags, Aerosol)
!  </TEMPLATE>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatology input
!  </IN>
!  <INOUT NAME="Aerosolrad_diags" TYPE="aerosolrad_diag_type">
!   Aerosol radiative diagnostics in radiation package
!  </INOUT>
!  <IN NAME="is, ie" TYPE="integer">
!   The longitude index of model physics window domain
!  </IN>
!  <IN NAME="js, je" TYPE="integer">
!   The latitude index of model physics window domain
!  </IN>
! </SUBROUTINE>
!
subroutine aerosol_radiative_properties (is, ie, js, je, &
                                         Time, fracday, p_half, relhum, deltaz, &
                                         do_cmip_sw_diagnostics, Aerosolrad_control, &
                                         Aerosolrad_diags, Aerosol, extinction, &
                                         aerooptdep, aerooptdep_volc, &
                                         aeroasymfac, aerosctopdep, aeroextopdep)

!---------------------------------------------------------------------
!    aerosol_radiative_properties defines and returns the radiative
!    properties for each aerosol properties type and for each solar
!    parameterization band in the shortwave and for each aerosol 
!    emissivity band in the longwave.
!---------------------------------------------------------------------

integer,                       intent(in)    :: is, ie, js, je
type(time_type),               intent(in)    :: Time
real, dimension(:,:),          intent(in)    :: fracday
real, dimension(:,:,:),        intent(in)    :: p_half, relhum, deltaz
logical,                       intent(in)    :: do_cmip_sw_diagnostics
type(aerosolrad_control_type), intent(in)    :: Aerosolrad_control
type(aerosol_type),            intent(in)    :: Aerosol
type(aerosolrad_diag_type),    intent(inout) :: Aerosolrad_diags
real, dimension(:,:,:),        intent(out)   :: extinction  !sw extinction for volcanoes
real, dimension(:,:,:,:),      intent(out)   :: aerooptdep, aerooptdep_volc, &
                                                aeroasymfac, aerosctopdep, aeroextopdep
 
!----------------------------------------------------------------------
! local variables:                                                     

      integer  :: na, nw, ni, n       ! do-loop indices
      integer  :: iaer, i, j, k
      integer  :: nextinct
      logical  :: including_swaerosols
      real, dimension (size(Aerosol%aerosol,1),  &
                       size(Aerosol%aerosol,2),  &
                       size(Aerosol%aerosol,3))  :: sul, bc
      integer, dimension (size(Aerosol%aerosol,1), &
                          size(Aerosol%aerosol,2), &
                          size(Aerosol%aerosol,3)) :: ivol

      logical, dimension (size(Aerosol%aerosol,1), &
                          size(Aerosol%aerosol,2)) :: daylight


!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosolrad_package_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate and initialize arrays to hold aerosol diagnostics.
!---------------------------------------------------------------------

      call Aerosolrad_diags%alloc ( ie-is+1, je-js+1, &
                                   size(Aerosol%aerosol,3), &
                                   size(Aerosol%aerosol,4), &
                                   nfields_sw_ext, &
                                   nfields_sw_ssa, &
                                   nfields_sw_asy, &
                                   nfields_lw_ext, &
                                   nfields_lw_ssa, &
                                   nfields_lw_asy, &
                                   Aerosolrad_control )
     !call aerosolrad_diag_alloc ( Aerosolrad_diags, &
     !                             ie-is+1, je-js+1, &
     !                             size(Aerosol%aerosol,3), &
     !                             size(Aerosol%aerosol,4), &
     !                             Aerosolrad_control)

!--------------------------------------------------------------------
!    if the volcanic sw aerosol extinction is being supplied, obtain
!    the appropriate data.
!--------------------------------------------------------------------
        if (using_sw_ext) then

!---------------------------------------------------------------------
!    if new sw extinction data is needed on this step, call interpolator
!    to obtain it.  if the data is not to be interpolated, save the
!    retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_sw_ext) then
            if (nfields_sw_ext >= 1) then
              call interpolator (Sw_aer_extopdep_interp, Volcano_Time, &
                                 p_half, Aerosolrad_diags%sw_ext,    &
                                 sw_ext_name(1), is, js)
            endif
            if (.not. interpolating_volcanic_data) then
              sw_ext_save(is:ie,js:je,:,:) = Aerosolrad_diags%sw_ext
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if ( .not. interpolating_volcanic_data) then
              Aerosolrad_diags%sw_ext = sw_ext_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic sw aerosol single scattering albedo is being 
!    supplied, obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_sw_ssa) then

!---------------------------------------------------------------------
!    if new sw single scattering albedo data is needed on this step, 
!    call interpolator to obtain it.  if the data is not to be inter-
!    polated, save the retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_sw_ssa) then
            if (nfields_sw_ssa >= 1) then
              call interpolator (Sw_aer_ssalb_interp, Volcano_Time, &
                                 p_half, Aerosolrad_diags%sw_ssa,    &
                                 sw_ssa_name(1), is, js)
            endif
            if ( .not. interpolating_volcanic_data) then
              sw_ssa_save(is:ie,js:je,:,:) = Aerosolrad_diags%sw_ssa
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if ( .not. interpolating_volcanic_data) then
              Aerosolrad_diags%sw_ssa = sw_ssa_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic sw aerosol asymmetry factor is being supplied, 
!    obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_sw_asy) then

!---------------------------------------------------------------------
!    if new sw asymmetry factor data is needed on this step, call 
!    interpolator to obtain it.  if the data is not to be interpolated,
!    save the retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_sw_asy) then
            if (nfields_sw_asy >= 1) then
              call interpolator (Sw_aer_asymm_interp, Volcano_Time, &
                                 p_half, Aerosolrad_diags%sw_asy,    &
                                 sw_asy_name(1), is, js)
            endif
            if ( .not. interpolating_volcanic_data) then
              sw_asy_save(is:ie,js:je,:,:) = Aerosolrad_diags%sw_asy
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if (.not. interpolating_volcanic_data) then
              Aerosolrad_diags%sw_asy = sw_asy_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw aerosol extinction is being supplied, obtain
!    the appropriate data.
!--------------------------------------------------------------------
        if (using_lw_ext) then

!---------------------------------------------------------------------
!    if new lw extinction data is needed on this step, call interpolator
!    to obtain it.  if the data is not to be interpolated, save the
!    retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_lw_ext) then
            if (nfields_lw_ext >= 1) then
              call interpolator (Lw_aer_extopdep_interp, Volcano_Time, &
                                 p_half, Aerosolrad_diags%lw_ext,    &
                                 lw_ext_name(1), is, js)
            endif
            if (.not. interpolating_volcanic_data) then
              lw_ext_save(is:ie,js:je,:,:) = Aerosolrad_diags%lw_ext
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if ( .not. interpolating_volcanic_data) then
              Aerosolrad_diags%lw_ext = lw_ext_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw single scattering albedo is being supplied, 
!    obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_lw_ssa) then

!---------------------------------------------------------------------
!    if new lw single scattering albedo data is needed on this step, 
!    call interpolator to obtain it.  if the data is not to be inter-
!    polated, save the retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_lw_ssa) then
            if (nfields_lw_ssa >= 1) then
              call interpolator (Lw_aer_ssalb_interp, Volcano_Time, &
                                 p_half, Aerosolrad_diags%lw_ssa,    &
                                 lw_ssa_name(1), is, js)
            endif
            if ( .not. interpolating_volcanic_data) then
              lw_ssa_save(is:ie,js:je,:,:) = Aerosolrad_diags%lw_ssa
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if ( .not. interpolating_volcanic_data) then
              Aerosolrad_diags%lw_ssa = lw_ssa_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw aerosol asymmetry factor is being supplied, 
!    obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_lw_asy) then

!---------------------------------------------------------------------
!    if new lw asymmetry factor data is needed on this step, call 
!    interpolator to obtain it.  if the data is not to be interpolated, 
!    save the retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_lw_asy) then
            if (nfields_lw_asy >= 1) then
              call interpolator (Lw_aer_asymm_interp, Volcano_Time, &
                                 p_half, Aerosolrad_diags%lw_asy,    &
                                 lw_asy_name(1), is, js)
            endif
            if (.not. interpolating_volcanic_data) then
              lw_asy_save(is:ie,js:je,:,:) = Aerosolrad_diags%lw_asy
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if (.not. interpolating_volcanic_data) then
              Aerosolrad_diags%lw_asy = lw_asy_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!---------------------------------------------------------------------
!    code for treating sulfate and black carbon as an internal aerosol
!    mixture.
!---------------------------------------------------------------------
        if (using_im_bcsul) then
          if (num_sul > 0) then
            sul(:,:,:) = Aerosol%aerosol(:,:,:,sul_ind(1))
            do iaer=2,num_sul
              sul(:,:,:) = sul(:,:,:) +    &
                                Aerosol%aerosol(:,:,:,sul_ind(iaer))
            end do
          else
            sul = 0.
          endif
          if (num_bc > 0) then
            bc(:,:,:) = Aerosol%aerosol(:,:,:,bc_ind(1))
            do iaer=2,num_bc
              bc(:,:,:) = bc(:,:,:) +    &
                                Aerosol%aerosol(:,:,:,bc_ind(iaer))
            end do
          else
            bc = 0.
          endif
          do k = 1,size(Aerosol%aerosol,3)
            do j = 1,size(Aerosol%aerosol,2)
              do i = 1,size(Aerosol%aerosol,1)
                if (bc(i,j,k) > 0 .and. sul(i,j,k) > 0.0) then
                  ivol(i,j,k) = 100-MIN(100, MAX( 0,     &
                   NINT(100.*sul(i,j,k)/(sul(i,j,k) +bc(i,j,k)*1.74))))
                else
                  ivol(i,j,k) = 0
                end if
              enddo
            end do
          end do
        else
          ivol = 0
        endif ! (using_im_bcsul)

!--------------------------------------------------------------------
!    longwave aerosol optical depth
!--------------------------------------------------------------------

      do n=1,N_AEROSOL_BANDS  !  loop on aerosol frequency bands
        if (Aerosolrad_control%do_lwaerosol .or. Aerosolrad_control%do_lwaerosol_forcing) then
          call lw_optical_depth_aerosol ( relhum, n, &
                                          Aerosol, ivol, Aerosolrad_diags, &
                                          aerooptdep(:,:,:,n) )
        else
          aerooptdep(:,:,:,n) = 0.0
        endif

        if (Aerosolrad_control%volcanic_lw_aerosols) then
          call lw_volcanic_aerosol ( deltaz, n, Aerosolrad_diags, &
                                     aerooptdep_volc(:,:,:,n) )
        else
          aerooptdep_volc(:,:,:,n) = 0.0
        endif
      end do
!--------------------------------------------------------------------
!    define a flag indicating columns in which there is sunshine during
!    this radiation time step. define a flag indicating points with both
!    sunlight and cloud.      
!----------------------------------------------------------------------c
      do j = 1, size(Aerosol%aerosol,2)
        do i = 1, size(Aerosol%aerosol,1)
          if ( fracday(i,j) /= 0.0) then 
            daylight(i,j) = .true.     
          else 
            daylight(i,j) = .false.     
          endif     
        end do
      end do

     ! calculate_volcanic_sw_heating in shortwave_driver_nml

      if (Aerosolrad_control%do_swaerosol .or. Aerosolrad_control%do_swaerosol_forcing) then
        including_swaerosols = .true.
      else
        including_swaerosols = .false.
      endif

     !if (calculate_volcanic_sw_heating) then
     !    call sw_aerosol_optical_props     &    
     !              (Atmos_input%aerosolrelhum, Atmos_input%deltaz,
     !               Aerosol, ivol, &
     !               Aerosolrad_control%volcanic_sw_aerosols, Aerosolrad_diags,  &
     !               Aerosolrad_control%do_swaerosol, do_cmip_sw_diagnostics, &
     !               daylight, aeroextopdep, aerosctopdep, aeroasymfac, &
     !               aeroextopdep2, aerosctopdep2, aeroasymfac2) 
     !else
          call sw_aerosol_optical_props     &    
                    (relhum, deltaz, &
                     Aerosol, ivol, &
                     Aerosolrad_control%volcanic_sw_aerosols, Aerosolrad_diags,  &
                     including_swaerosols, do_cmip_sw_diagnostics, &
                !BW  .true., do_cmip_sw_diagnostics, &
                     daylight, aeroextopdep, aerosctopdep, aeroasymfac)
     !endif

     ! SW extinction (band 4) for vocanoes
     ! this will be used in the atmos tracer driver
     extinction(:,:,:) = 0.0
     if (Aerosolrad_control%volcanic_sw_aerosols .and. Aerosolrad_control%do_swaerosol) then 
         nextinct = get_tracer_index(MODEL_ATMOS,'Extinction')
         if (nextinct /= NO_TRACER) extinction(:,:,:) = Aerosolrad_diags%sw_ext(:,:,:,4)
     endif

!--------------------------------------------------------------------

end subroutine aerosol_radiative_properties

!#####################################################################
! <SUBROUTINE NAME="aerosolrad_package_end">
!  <OVERVIEW>
!    aerosolrad_package_end is the destructor for 
!    aerosolrad_package_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    aerosolrad_package_end is the destructor for 
!    aerosolrad_package_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosolrad_package_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine aerosolrad_package_end (Aerosolrad_control)

type(aerosolrad_control_type), intent(in) :: Aerosolrad_control

!--------------------------------------------------------------------
!    aerosolrad_package_end is the destructor for 
!    aerosolrad_package_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosolrad_package_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    deallocate module variables.
!---------------------------------------------------------------------
      if (Aerosolrad_control%do_swaerosol .or. Aerosolrad_control%do_swaerosol_forcing) then      
        deallocate (solivlaero, nivl1aero, nivl2aero, endaerwvnsf, &
                    aeroextivl, aerossalbivl, aeroasymmivl)
      endif
      if (Aerosolrad_control%do_lwaerosol .or. Aerosolrad_control%do_lwaerosol_forcing) then
        deallocate ( sflwwts)
        deallocate ( sflwwts_cn)
      endif
      
!---------------------------------------------------------------------
!    deallocate module arrays
!---------------------------------------------------------------------
      if (allocated(sulfate_index_MOD )) deallocate(sulfate_index_MOD )
      if (allocated(bcphilic_index_MOD)) deallocate(bcphilic_index_MOD)
      if (allocated(omphilic_index_MOD)) deallocate(omphilic_index_MOD)
      if (allocated(seasalt1_index_MOD)) deallocate(seasalt1_index_MOD)
      if (allocated(seasalt2_index_MOD)) deallocate(seasalt2_index_MOD)
      if (allocated(seasalt3_index_MOD)) deallocate(seasalt3_index_MOD)
      if (allocated(seasalt4_index_MOD)) deallocate(seasalt4_index_MOD)
      if (allocated(seasalt5_index_MOD)) deallocate(seasalt5_index_MOD)
      if (allocated(seasalt_aitken_index_MOD)) deallocate(seasalt_aitken_index_MOD)
      if (allocated(seasalt_fine_index_MOD)) deallocate(seasalt_fine_index_MOD)
      if (allocated(seasalt_coarse_index_MOD)) deallocate(seasalt_coarse_index_MOD)
      if (allocated(nitrate_index_MOD)) deallocate(nitrate_index_MOD)
      if (allocated(optical_index_MOD )) deallocate(optical_index_MOD )

      
      if (Aerosolrad_control%volcanic_lw_aerosols) then
        if (nfields_lw_ext /= 0) then
          call interpolator_end (Lw_aer_extopdep_interp)
        endif
        if (nfields_lw_ssa /= 0) then
        call interpolator_end (Lw_aer_ssalb_interp)
        endif
        if (nfields_lw_asy /= 0) then
        call interpolator_end (Lw_aer_asymm_interp)
        endif
      endif

      if (Aerosolrad_control%volcanic_sw_aerosols) then
        if (nfields_sw_ext /= 0) then
        call interpolator_end (Sw_aer_extopdep_interp)
        endif
        if (nfields_sw_ssa /= 0) then
        call interpolator_end (Sw_aer_ssalb_interp)
        endif
        if (nfields_sw_asy /= 0) then
        call interpolator_end (Sw_aer_asymm_interp)
        endif
      endif

      if (.not. interpolating_volcanic_data) then
        deallocate (sw_ext_save, sw_ssa_save, sw_asy_save, &
                    lw_ext_save, lw_ssa_save, lw_asy_save)
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

end subroutine aerosolrad_package_end

!####################################################################

 subroutine number_of_lw_aerosol_bands (nbands)

!-----------------------------------------------------------------------
!    number_of_lw_aerosol_bands returns the total number of
!    IR aerosol emissivity bands.
!-----------------------------------------------------------------------

 integer,  intent(out) :: nbands

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'aerosolrad_package_mod',  &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------

      nbands = N_AEROSOL_BANDS

!---------------------------------------------------------------------

 end subroutine number_of_lw_aerosol_bands



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  
                                  

!#####################################################################
! <SUBROUTINE NAME="assign_aerosol_opt_props">
!  <OVERVIEW>
!    assign_aerosol_opt_props assigns an index for an available optical
!    properties type to each activated aerosol type. for sulfates, a 
!    flag is set, since the aerosol properties type is a function 
!    of model relative humidity, and will vary with time.
!  </OVERVIEW>
!  <DESCRIPTION>
!    assign_aerosol_opt_props assigns an index for an available optical
!    properties type to each activated aerosol type. for sulfates, a 
!    flag is set, since the aerosol properties type is a function 
!    of model relative humidity, and will vary with time.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call assign_aerosol_opt_props (aerosol_names)
!  </TEMPLATE>
!  <IN NAME="aerosol_names" TYPE="character">
!   names associated with each aerosol species
!  </IN>
! </SUBROUTINE>
!
subroutine assign_aerosol_opt_props (aerosol_names)

!----------------------------------------------------------------------
!    assign_aerosol_opt_props assigns an index for an available optical
!    properties type to each activated aerosol type. for sulfates, a 
!    flag is set, since the aerosol properties type is a function 
!    of model relative humidity, and will vary with time.
!---------------------------------------------------------------------

!character(len=64), dimension(:), intent(in) :: aerosol_names
character(len=*), dimension(:), intent(in) :: aerosol_names

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     aerosol_names     names associated with each aerosol species
!
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    local variables:

      character(len=64) :: name_in, target_name
!yim
      character(len=4)  :: chind, chind2
      integer           :: nfields
!yim
      integer           :: n, noptical, m
      integer           :: ibc, isul

!---------------------------------------------------------------------
!   local variables:
!
!       name_in          variable to hold current aerosol name 
!                        being processed
!       target_name      aerosol_optical_name associated with a given
!                        aerosol species     
!       nfields          number of activated aerosol species
!       n, noptical      do-loop indices
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    count the number of aerosol optical property categories requested
!    via the namelist input.
!----------------------------------------------------------------------
      do n=1,MAX_OPTICAL_FIELDS
        if (aerosol_optical_names(n) /= ' '  ) then
          naermodels = n
        else
          exit
        endif
      end do !n

!---------------------------------------------------------------------
!    define the number of activated aerosol species.
!---------------------------------------------------------------------
      nfields = size (aerosol_names(:))

      if (using_im_bcsul) then
        allocate (sul_ind(nfields))
        allocate (bc_ind(nfields))
      endif

!---------------------------------------------------------------------
!    allocate module arrays which will contain the indices
!    for the different aerosol mixtures
!---------------------------------------------------------------------
      if (using_im_bcsul) then
        allocate (sulfate_index_MOD (0:100,0:100))
      else
        allocate (sulfate_index_MOD (0:100,0:0  ))
      endif
      allocate (omphilic_index_MOD(0:100 ), &
                bcphilic_index_MOD(0:100 ), &
                seasalt1_index_MOD(0:100 ), &
                seasalt2_index_MOD(0:100 ), &
                seasalt3_index_MOD(0:100 ), &
                seasalt4_index_MOD(0:100 ), &
                seasalt5_index_MOD(0:100 ), & 
                seasalt_aitken_index_MOD(0:100 ), &
                seasalt_fine_index_MOD(0:100 ), &
                seasalt_coarse_index_MOD(0:100 ), &
                nitrate_index_MOD(0:100 ) )
      allocate (optical_index_MOD(nfields) )
      optical_index_MOD    = 0
      sulfate_index_MOD    = 0
      omphilic_index_MOD    = 0
      bcphilic_index_MOD    = 0
      seasalt1_index_MOD    = 0
      seasalt2_index_MOD    = 0
      seasalt3_index_MOD    = 0
      seasalt4_index_MOD    = 0
      seasalt5_index_MOD    = 0
      seasalt_aitken_index_MOD = 0
      seasalt_fine_index_MOD   = 0
      seasalt_coarse_index_MOD = 0
      nitrate_index_MOD    = 0

!----------------------------------------------------------------------
!    match aerosol optical property indices with aerosol indices.
!    sulfate aerosols are handled separately (below) with RH dependence.
!----------------------------------------------------------------------
      num_sul = 0
      num_bc = 0
      isul = 1
      ibc = 1
      do n=1,nfields
        name_in = trim(aerosol_names(n))

        if (name_in == 'so4' .or. name_in == 'so4_anthro' .or. name_in == 'so4_natural') then
          optical_index_MOD(n) = SULFATE_FLAG
          if (using_im_bcsul) then
            num_sul = num_sul +1
            sul_ind(isul) = n
            isul = isul + 1
          endif
        else if (name_in == "omphilic" .or. name_in == "oc_hydrophilic") then
           optical_index_MOD(n) = OMPHILIC_FLAG
        else if (name_in == "bcphilic" .or. name_in == "bc_hydrophilic") then
           optical_index_MOD(n) = BCPHILIC_FLAG
           if (using_im_bcsul) then
              num_bc = num_bc +1
              bc_ind(ibc) = n
              ibc = ibc + 1
           endif
        else if (name_in == "nitrate") then
            optical_index_MOD(n) = NITRATE_FLAG
        else if (name_in == "seasalt1") then
            optical_index_MOD(n) = SEASALT1_FLAG
        else if (name_in == "seasalt2") then
            optical_index_MOD(n) = SEASALT2_FLAG
        else if (name_in == "seasalt3") then
            optical_index_MOD(n) = SEASALT3_FLAG
        else if (name_in == "seasalt4") then
            optical_index_MOD(n) = SEASALT4_FLAG
        else if (name_in == "seasalt5") then
            optical_index_MOD(n) = SEASALT5_FLAG
        else if (name_in == "seasalt_aitken") then
            optical_index_MOD(n) = SEASALTA_FLAG
        else if (name_in == "seasalt_fine") then
            optical_index_MOD(n) = SEASALTF_FLAG
        else if (name_in == "seasalt_coarse") then
            optical_index_MOD(n) = SEASALTC_FLAG
!yim
        else if (name_in == "black_carbon" .and. using_im_bcsul) then
            optical_index_MOD(n) = BC_FLAG
            num_bc = num_bc +1
            bc_ind(ibc) = n
            ibc = ibc + 1
        else 
          select case( name_in )
            case( "anthro_dust_0.1", "natural_dust_0.1" )
              target_name = "dust_0.1"
            case( "anthro_dust_0.2", "natural_dust_0.2" )
              target_name = "dust_0.2"
            case( "anthro_dust_0.4", "natural_dust_0.4" )
              target_name = "dust_0.4"
            case( "anthro_dust_0.8", "natural_dust_0.8" )
              target_name = "dust_0.8"
            case( "anthro_dust_1.0", "natural_dust_1.0" )
              target_name = "dust_1.0"
            case( "anthro_dust_2.0", "natural_dust_2.0" )
              target_name = "dust_2.0"
            case( "anthro_dust_4.0", "natural_dust_4.0" )
              target_name = "dust_4.0"
            case( "anthro_dust_8.0", "natural_dust_8.0" )
              target_name = "dust_8.0"
            case( "black_carbon" )
               target_name = "soot"
            case( "organic_carbon" )
              target_name = "organic_carbon"
            case( "sea_salt" )
              target_name = "sea_salt"
            case( "dust1" )
              target_name = "dust1"
            case( "dust2" )
              target_name = "dust2"
            case( "dust3" )
              target_name = "dust3"
            case( "dust4" )
              target_name = "dust4"
            case( "dust5" )
              target_name = "dust5"
            case( "bcdry" )
              target_name = "bcdry"
            case( "omphobic", "oc_hydrophobic" )
              target_name = "omphobic"
            case( "bcphobic", "bc_hydrophobic" )
              target_name = "bcphobic"
            case( "dust_mode1_of_1" )
              target_name = "dust_mode1_of_1"
            case( "dust_mode1_of_2" )
              target_name = "dust_mode1_of_2"
            case( "dust_mode2_of_2" )
              target_name = "dust_mode2_of_2"
            case( "dust_mode1_of_3" )
              target_name = "dust_mode1_of_3"
            case( "dust_mode2_of_3" )
              target_name = "dust_mode2_of_3"
            case( "dust_mode3_of_3" )
              target_name = "dust_mode3_of_3"
            case DEFAULT
              target_name = name_in
          end select  

!--------------------------------------------------------------------
!    go through the set of aerosol properties types looking for 
!    the target_name defined above. when found, associate the
!    optical properties type index with the current aerosol species.
!--------------------------------------------------------------------
          do noptical=1,naermodels
            if (aerosol_optical_names(noptical) == target_name) then
              optical_index_MOD(n) = noptical
              exit
            end if
          end do !noptical

!--------------------------------------------------------------------
!    if the target_name is not found, exit with an error message.
!----------------------------------------------------------------------
          if (optical_index_MOD(n) == 0 ) then
            call error_mesg( 'aerosolrad_package_mod', &
                'Cannot find aerosol optical model = ' //    &
                                           TRIM( target_name ), FATAL )
          endif
        endif  ! (name_in ==)
      end do  ! (n=1,nfields)

  select case(trim(aerosol_data_set))
    case ('Ginoux_Reddy') 

     if (using_im_bcsul) then

!----------------------------------------------------------------------
!    set up RH-dependent sulfate aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
       do n=0,100
          do m=0,100
            chind = integer_string(sulfate_indices(n))
            chind2 = integer_string(sulfate_vol_indices(m))
!yim format sulfate_10%_10% (RH + volume fraction)
            target_name = 'sulfate_' // trim(chind)  // '%_' // trim(chind2)// '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
            do noptical=1,naermodels
               if (aerosol_optical_names(noptical) == target_name ) then
                  sulfate_index_MOD(n,m) = noptical
                  exit 
               end if
            end do ! noptical

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
            if (sulfate_index_MOD(n,m) == 0 ) then
               call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
            endif
        end do !m
      end do !n
     else ! (.not.) using_im_bcsul
!----------------------------------------------------------------------
!    set up RH-dependent sulfate aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
       call optical_property_index ('sulfate', sulfate_indices, sulfate_index_MOD(:,0))

     endif !using_im_bssul

!---------------------------------------------------------------------
!    set up RH-dependent nitrate aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
     call optical_property_index ('nitrate', nitrate_indices, nitrate_index_MOD)

!---------------------------------------------------------------------
!    set up RH-dependent omphilic aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
     call optical_property_index ('omphilic', omphilic_indices, omphilic_index_MOD)

!---------------------------------------------------------------------
!    set up RH-dependent seasalt1 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
     call optical_property_index ('seasalt1', seasalt1_indices, seasalt1_index_MOD)

!---------------------------------------------------------------------
!    set up RH-dependent seasalt2 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
     call optical_property_index ('seasalt2', seasalt2_indices, seasalt2_index_MOD)

!---------------------------------------------------------------------
!    set up RH-dependent seasalt3 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
     call optical_property_index ('seasalt3', seasalt3_indices, seasalt3_index_MOD)

!---------------------------------------------------------------------
!    set up RH-dependent seasalt4 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
     call optical_property_index ('seasalt4', seasalt4_indices, seasalt4_index_MOD)

!-------------------------------------------------------------------
!    set up RH-dependent seasalt5 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
     call optical_property_index ('seasalt5', seasalt5_indices, seasalt5_index_MOD)

!---------------------------------------------------------------------
!    set up RH-dependent seasalt_aitken aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------

     call optical_property_index('seasalt_aitken',seasalt_aitken_indices, seasalt_aitken_index_MOD)

!---------------------------------------------------------------------
!    set up RH-dependent seasalt_fine aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------

     call optical_property_index('seasalt_fine', seasalt_fine_indices, seasalt_fine_index_MOD)

!---------------------------------------------------------------------
!    set up RH-dependent seasalt_coarse aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------

     call optical_property_index('seasalt_coarse', seasalt_coarse_indices, seasalt_coarse_index_MOD)

!-------------------------------------------------------------------
!    set up RH-dependent bcphilic aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
     if ( .not. using_im_bcsul )  then
          call optical_property_index ('bcphilic', bcphilic_indices, bcphilic_index_MOD)
     endif


  case ('shettle_fenn')

   if (using_im_bcsul) then
!----------------------------------------------------------------------
!    set up RH-dependent sulfate aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
          do m=0,100
            chind = integer_string(sulfate_indices(n))
            chind2 = integer_string(sulfate_vol_indices(m))
!yim format sulfate_10%_10% (RH + volume fraction)
            target_name = 'sulfate_' // trim(chind)  // '%_' // trim(chind2)// '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == target_name ) then
            sulfate_index_MOD(n,m) = noptical
            exit
          end if
        end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
        if (sulfate_index_MOD(n,m) == 0 ) then
          call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
        endif
      end do
      end do
   else 
!----------------------------------------------------------------------
!    set up RH-dependent sulfate aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          call optical_property_index ('sulfate', sulfate_indices, sulfate_index_MOD(:,0))

   endif
  end select  ! (aerosol_data_set)  

!---------------------------------------------------------------------


end subroutine assign_aerosol_opt_props


!######################################################################
! <SUBROUTINE NAME="read_optical_input_file">
!  <OVERVIEW>
!    read_optical_input_file reads the optical properties input file
!    to obtain the specified aerosol radiative properties for each 
!    aerosol in each of the aerosol parameterization spectral bands.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_optical_input_file reads the optical properties input file
!    to obtain the specified aerosol radiative properties for each 
!    aerosol in each of the aerosol parameterization spectral bands.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_optical_input_file
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine read_optical_input_file

!-----------------------------------------------------------------------
!    read_optical_input_file reads the optical properties input file
!    to obtain the specified aerosol radiative properties for each 
!    aerosol in each of the aerosol parameterization spectral bands.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    local variables:

      real,    dimension(:), allocatable    :: aeroext_in,   &
                                               aerossalb_in,   &
                                               aeroasymm_in
      logical, dimension(:), allocatable    :: found

      integer           :: unit, num_input_categories
      character(len=64) :: name_in
      integer           :: n, noptical

!---------------------------------------------------------------------
!   local variables:
!
!       aeroext_in       aerosol extinction coefficient read from 
!                        input file
!       aerossalb_in     aerosol single scattering albedo read from 
!                        input file
!       aeroasymm_in     aerosol asymmetry factor read from 
!                        input file
!       found            aerosol radiative property data has been
!                        obtained from input file for the given
!                        optical properties type ?
!       unit             io unit number used for optical properties file
!       num_input_categories
!                        number of optical properties types contained
!                        in optical data input file
!       name_in          name of optical properties type being processed
!       n, noptical      do-loop indices
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    open the ASCII input file containing aerosol optical property
!    information.
!----------------------------------------------------------------------
      call mpp_open (unit, 'INPUT/'//optical_filename, MPP_RDONLY,  &
                     MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)

!----------------------------------------------------------------------
!    read the dimension information contained in the input file.
!----------------------------------------------------------------------
      read ( unit,* ) num_wavenumbers
      read ( unit,* ) num_input_categories

!----------------------------------------------------------------------
!    read wavenumber limits for aerosol parameterization bands from 
!    the input file.
!----------------------------------------------------------------------
       allocate (endaerwvnsf(num_wavenumbers) )
       read (unit,* )
       read (unit,* ) endaerwvnsf
 
!----------------------------------------------------------------------
!    allocate module arrays to hold the specified sw properties for 
!    each parameterization bnad and each aerosol properties type.
!----------------------------------------------------------------------
      allocate (       &
            aeroextivl   (num_wavenumbers, naermodels),&
            aerossalbivl (num_wavenumbers, naermodels), &
            aeroasymmivl (num_wavenumbers, naermodels) )

!----------------------------------------------------------------------
!    allocate local working arrays.
!----------------------------------------------------------------------
      allocate (aeroext_in   (num_wavenumbers ),             &
                aerossalb_in (num_wavenumbers ),           &
                aeroasymm_in (num_wavenumbers ),           &
                found        (naermodels ) )

!----------------------------------------------------------------------
!    match the names of optical property categories from input file with
!    those specified in the namelist, and store the following data
!    appropriately. indicate that the data has been found.
!----------------------------------------------------------------------
      found(:) = .false.
      do n=1,num_input_categories
        read( unit,* ) name_in
        read( unit,* )
        read( unit,* ) aeroext_in
        read( unit,* )
        read( unit,* ) aerossalb_in
        read( unit,* )
        read( unit,* ) aeroasymm_in
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == name_in) then
            aeroextivl(:,noptical)   = aeroext_in
            aerossalbivl(:,noptical) = aerossalb_in
            aeroasymmivl(:,noptical) = aeroasymm_in
            found( noptical ) = .true.
            exit
          endif
        end do
      end do

!----------------------------------------------------------------------
!    close the ASCII input file.
!----------------------------------------------------------------------
      call mpp_close( unit )

!----------------------------------------------------------------------
!    check to make sure data for all aerosol optical property
!    categories specified in namelist were contained in ASCII
!    input file. if not, exit with a message.
!----------------------------------------------------------------------
      do noptical = 1,naermodels
        if (.not. found( noptical ) ) then
              call error_mesg( 'aerosolrad_package_mod', &
              'Cannot find aerosol optical properties for ' // &
                TRIM(aerosol_optical_names(noptical)),  FATAL )
        endif
      end do

!----------------------------------------------------------------------
!    deallocate local working arrays.
!----------------------------------------------------------------------
      deallocate (aeroext_in, aerossalb_in, aeroasymm_in, found)



end subroutine read_optical_input_file


!#####################################################################

function integer_string (ivalue)
integer, intent(in) :: ivalue
character(len=4) :: integer_string

!-----------------------------------------------------------------------
!    integer_string encodes an integer value into an integer string
!    integer values can range from 0 to 100.
!-----------------------------------------------------------------------
      if (ivalue >= 0 .and. ivalue < 10) then
        write (integer_string, '(i1)') ivalue
      else if (ivalue >= 10 .and. ivalue < 100) then
        write (integer_string, '(i2)') ivalue
      else if (ivalue == 100) then
        write (integer_string, '(i3)') ivalue
      else
        call error_mesg ('integer_string', &
                   'Integer value must be between 0 and 100', FATAL)
      endif

!-----------------------------------------------------------------------

end function integer_string

!#####################################################################

subroutine optical_property_index (aerosol_name, aerosol_indices, &
                                   aerosol_index)

character(len=*), intent(in)                      :: aerosol_name
integer,          intent(in),  dimension(0:100)   :: aerosol_indices
integer,          intent(out), dimension(0:100)   :: aerosol_index

!---------------------------------------------------------------------
integer           ::  n, noptical
character(len=4)  :: chind
character(len=64) :: target_name

!---------------------------------------------------------------------
!    set up RH-dependent aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
      do n=0,100
        chind = integer_string(aerosol_indices(n))
        target_name = trim(aerosol_name) // '_' // trim(chind)  // '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == target_name ) then
            aerosol_index(n) = noptical
            exit
          endif
        end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
        if (aerosol_index(n) == 0 ) then
          call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                        TRIM( target_name), FATAL )
        endif

      end do
!---------------------------------------------------------------------

end subroutine optical_property_index

!#####################################################################
! <SUBROUTINE NAME="sw_aerosol_interaction">
!  <OVERVIEW>
!    sw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the solar parameterization
!    bands being used by the model.
!  </OVERVIEW>
!  <DESCRIPTION>
!    sw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the solar parameterization
!    bands being used by the model.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call sw_aerosol_interaction
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine sw_aerosol_interaction

!-----------------------------------------------------------------------
!    sw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the solar parameterization
!    bands being used by the model.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!    local variables:

      integer           :: nbands, nband, nivl3
      real              :: sumsol3
      integer           :: nw
      integer           :: nmodel

!---------------------------------------------------------------------
!   local variables:
!
!       aeroext_in       aerosol extinction coefficient read from 
!                        input file
!       aerossalb_in     aerosol single scattering albedo read from 
!                        input file
!       aeroasymm_in     aerosol asymmetry factor read from 
!                        input file
!       found            aerosol radiative property data has been
!                        obtained from input file for the given
!                        optical properties type ?
!       unit             io unit number used for optical properties file
!       num_input_categories
!                        number of optical properties types contained
!                        in optical data input file
!       name_in          name of optical properties type being processed
!       nbands           number of bands in solar spectral param-
!                        eterization
!       nband            currently active solar spectrum band 
!       nivl3            currently active aerosol parameterization band
!       sumsol3          sum of solar input in current aerosol param-
!                        eterization band
!       n, nw, noptical  do-loop indices
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!    define the number of bands in the solar spectrum parameterization.
!    allocate space for variables defining the highest and lowest 
!    aerosol parameterization wavenumber in each solar spectral band
!    and the solar flux common to solar spectral band n and aerosol
!    parameterization band ni.
!---------------------------------------------------------------------
      call esfsw_number_of_bands (nbands)
      allocate ( nivl1aero  (nbands) )
      allocate ( nivl2aero  (nbands) )
      allocate ( solivlaero (nbands, num_wavenumbers))

!---------------------------------------------------------------------
!    define the solar weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the solar
!    spectral intervals and so determine the single-scattering proper-
!    ties on the solar spectral intervals.
!--------------------------------------------------------------------
      call esfsw_band_segments (endaerwvnsf, solivlaero, nivl1aero, nivl2aero)

!---------------------------------------------------------------------
!    allocate and initialize variables which will hold the aerosol 
!    radiative properties for each solar spectral parameterization band.
!    aerextband     the solar band values of the extinction 
!                   coefficient for aerosols                           
!    aerssalbband   the solar band values of the single-     
!                   scattering albedo for aerosols                      
!    aerasymmband   the solar band values of the asymmetry   
!                   factor for aerosols                                 
!---------------------------------------------------------------------
      allocate     &
        (aerextband_MOD   (nbands, naermodels), &
         aerssalbband_MOD (nbands, naermodels), &
         aerasymmband_MOD (nbands, naermodels) )
      aerextband_MOD   = 0.
      aerssalbband_MOD = 0.
      aerasymmband_MOD = 0.

!--------------------------------------------------------------------
!    if sw aerosol properties are desired and have not yet been calc-
!    ulated, use the thick-averaging technique to define the single-
!    scattering properties for each solar parameterization band n 
!    from the specified properties on the aerosol parameterization 
!    bands ni for each aerosol properties type nmodel. 
! references:                                                          
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!    code I: choosing a configuration for a large-scale model.,     
!    q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: a thin-averaging technique (subroutine esfsw_thinavg in 
!    esfsw_bands_mod) is also available.   
!--------------------------------------------------------------------
      do nmodel=1,naermodels
        call esfsw_thickavg (nivl1aero, nivl2aero, num_wavenumbers,   &
                       nbands, aeroextivl(:,nmodel), &
                       aerossalbivl(:,nmodel),    &
                       aeroasymmivl(:,nmodel), solivlaero,   &
                       aerextband_MOD(:,nmodel),    &
                       aerssalbband_MOD(:,nmodel),   &
                       aerasymmband_MOD(:,nmodel))
      end do

!---------------------------------------------------------------------


end subroutine sw_aerosol_interaction   


!#####################################################################
! <SUBROUTINE NAME="lw_aerosol_interaction">
!  <OVERVIEW>
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!  </OVERVIEW>
!  <DESCRIPTION>
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_aerosol_interaction
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine lw_aerosol_interaction      

!----------------------------------------------------------------------
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables:

!---------------------------------------------------------------------
!    the following arrays define the wavenumber ranges for the separate
!    aerosol emissivity bands in the model infrared parameterization. 
!    these may be changed only by the keeper of the radiation code.
!    the order of the frequency bands corresponds to the order used
!    in the lw radiation code.
!
!      aerbandlo_fr      low wavenumber limit for the non-continuum 
!                        aerosol emissivity bands
!      aerbandhi_fr      high wavenumber limit for the non-continuum
!                        aerosol emissivity bands
!      istartaerband_fr  starting wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      iendaerband_fr    ending wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      aerbandlo_co      low wavenumber limit for the continuum 
!                        aerosol emissivity bands
!      aerbandhi_co      high wavenumber limit for the continuum
!                        aerosol emissivity bands
!      istartaerband_co  starting wavenumber index for the continuum
!                        aerosol emissivity bands
!      iendaerband_co    ending wavenumber index for the continuum
!                        aerosol emissivity bands
!      aerbandlo         low wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      aerbandhi         high wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      istartaerband     starting wavenumber index for the entire set of
!                        aerosol emissivity bands
!      iendaerband       ending wavenumber index for the entire set of
!                        aerosol emissivity bands
!
!----------------------------------------------------------------------
      real, dimension (N_AEROSOL_BANDS_FR)     :: aerbandlo_fr =  &
      (/ 560.0, 630.0, 700.0, 800.0, 900.0,  990.0, 1070.0, 1200.0 /)

      real, dimension (N_AEROSOL_BANDS_FR)     :: aerbandhi_fr =  &
      (/ 630.0, 700.0, 800.0, 900.0, 990.0, 1070.0, 1200.0, 1400.0 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: istartaerband_fr =  &
      (/ 57,  64,  71,  81,  91, 100, 108, 121 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: iendaerband_fr =  &
      (/ 63,  70,  80,  90,  99, 107, 120, 140 /)

      real, dimension (N_AEROSOL_BANDS_CO)     :: aerbandlo_co =  &
      (/ 560.0 /)

      real, dimension (N_AEROSOL_BANDS_CO)     :: aerbandhi_co =  &
      (/ 800.0 /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: istartaerband_co =  &
      (/ 57  /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: iendaerband_co =  &
      (/ 80  /)

      integer, dimension (N_AEROSOL_BANDS_CN)  :: istartaerband_cn =  &
      (/ 81  /)

      integer, dimension (N_AEROSOL_BANDS_CN)  :: iendaerband_cn =  &
      (/ 120 /)

      real,    dimension(N_AEROSOL_BANDS)      :: aerbandlo, aerbandhi
      integer, dimension(N_AEROSOL_BANDS)      :: istartaerband,    &
                                                  iendaerband

!---------------------------------------------------------------------
!    the following arrays define how the ir aerosol band structure 
!    relates to the aerosol parameterization bands.
!
!      nivl1aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl2aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl1aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl2aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl1aer(n)       aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber for the 
!                        ir aerosol emissivity band n
!      nivl2aer(n)       aerosol parameterization band index corres-
!                        ponding to the highest wavenumber for the 
!                        ir aerosol emissivity band n
!      planckaerband(n)  planck function summed over each lw param-
!                        eterization band that is contained in the 
!                        ir aerosol emissivity band n
!
!---------------------------------------------------------------------
      integer, dimension (N_AEROSOL_BANDS_FR)  :: nivl1aer_fr,   &
                                                  nivl2aer_fr
      integer, dimension (N_AEROSOL_BANDS_CO)  :: nivl1aer_co,   &
                                                  nivl2aer_co
      integer, dimension (N_AEROSOL_BANDS_CN)  :: nivl1aer_cn,   &
                                                  nivl2aer_cn
      real,    dimension (N_AEROSOL_BANDS)     :: planckaerband
      real,    dimension (N_AEROSOL_BANDS_CN)  :: planckaerband_cn

!----------------------------------------------------------------------
!    the following arrays relate the ir aerosol emissivity band n to
!    either the aerosol optical properties type na or to the aerosol 
!    parameterization band ni.
!        aerextbandlw_fr(n,na)  band averaged extinction coefficient
!                               for non-continuum aerosol emissivity 
!                               band n and aerosol properties type na
!        aerssalbbandlw_fr(n,na)
!                               band averaged single-scattering
!                               coefficient for non-continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        aerextbandlw_co(n,na)  band averaged extinction coefficient
!                               for the continuum aerosol emissivity
!                               band n and aerosol properties type na
!        aerssalbbandlw_co(n,na)
!                               band averaged single-scattering
!                               coefficient for continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        planckivlaer_fr(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity non-
!                               continuum band n and aerosol parameter-
!                               ization band ni
!        planckivlaer_co(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity continuum 
!                               band n and aerosol parameterization 
!                               band ni
!        sflwwts_fr(n,ni)       band weights for the aerosol emissivity
!                               non-continuum band n and the aerosol 
!                               parameterization band ni 
!        sflwwts_co(n,ni)       band weights for the aerosol emissivity
!                               continuum band n and the aerosol 
!                               parameterization band ni 
!        planckivlaer(n,ni)     planck function over the spectral range
!                               common to aerosol emissivity band n and
!                               aerosol parameterization band ni
!        iendsfbands(ni)        ending wavenumber index for aerosol 
!                               parameterization band ni
!
!----------------------------------------------------------------------
      real,    dimension (N_AEROSOL_BANDS_FR, num_wavenumbers) :: &
                                                  planckivlaer_fr, &
                                                  sflwwts_fr
      real,    dimension (N_AEROSOL_BANDS_CO, num_wavenumbers) :: &
                                                  planckivlaer_co, &
                                                  sflwwts_co
      real,    dimension (N_AEROSOL_BANDS_CN, num_wavenumbers) :: &
                                                  planckivlaer_cn   
      integer, dimension (num_wavenumbers)    ::  iendsfbands

!---------------------------------------------------------------------
!    variables associated with the planck function calculation.
!    the planck function is defined for each of the NBLW longwave 
!    parameterization bands.
!---------------------------------------------------------------------
      real, dimension(NBLW)  :: c1, centnb, sc, src1nb, x, x1
      real                   :: del, xtemv, sumplanck

!---------------------------------------------------------------------
!    miscellaneous variables:

     logical         :: do_band1   !  should we do special calculation 
                                   !  for band 1 ?
     integer         :: ib, nw, nivl, nband, n, ni 
                                   !  do-loop indices and counters
     integer         :: na

!--------------------------------------------------------------------
!    define arrays containing the characteristics of all the ir aerosol
!    emissivity bands, both continuum and non-continuum.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        aerbandlo(n)     = aerbandlo_fr(n)
        aerbandhi(n)     = aerbandhi_fr(n)
        istartaerband(n) = istartaerband_fr(n)
        iendaerband(n)   = iendaerband_fr(n)
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        aerbandlo(n)     = aerbandlo_co     (n - N_AEROSOL_BANDS_FR)
        aerbandhi(n)     = aerbandhi_co     (n - N_AEROSOL_BANDS_FR)
        istartaerband(n) = istartaerband_co (n - N_AEROSOL_BANDS_FR)
        iendaerband(n)   = iendaerband_co   (n - N_AEROSOL_BANDS_FR)
      end do

!--------------------------------------------------------------------
!    allocate a module variable which will store the weighting function
!    between the aerosol emissivity bands and the aerosol parameter-
!    ization bands.
!--------------------------------------------------------------------
      allocate (sflwwts (N_AEROSOL_BANDS, num_wavenumbers))
      allocate (sflwwts_cn (N_AEROSOL_BANDS_CN, num_wavenumbers))

!--------------------------------------------------------------------
!    define the ending aerosol band index for each of the aerosol
!    parameterization bands.
!--------------------------------------------------------------------
      iendsfbands(:) = INT((endaerwvnsf(:) + 0.01)/10.0)

!--------------------------------------------------------------------
!    compute the planck function at 10C over each of the longwave
!    parameterization bands to be used as the weighting function. 
!--------------------------------------------------------------------
      do n=1,NBLW 
        del  = 10.0E+00
        xtemv = 283.15
        centnb(n) = 5.0 + (n - 1)*del
        c1(n)     = (3.7412E-05)*centnb(n)**3
        x(n)      = 1.4387E+00*centnb(n)/xtemv
        x1(n)     = EXP(x(n))
        sc(n)     = c1(n)/(x1(n) - 1.0E+00)
        src1nb(n) = del*sc(n)
      end do
 
!--------------------------------------------------------------------
!    sum the weighting function calculated over the longwave param-
!    eterization bands that are contained in each of the aerosol 
!    emissivity bands. 
!--------------------------------------------------------------------
      planckaerband(:) = 0.0E+00
      do n = 1,N_AEROSOL_BANDS
        do ib = istartaerband(n),iendaerband(n)
          planckaerband(n) = planckaerband(n) + src1nb(ib)
        end do
      end do
      planckaerband_cn(:) = 0.0E+00
      do n = 1,N_AEROSOL_BANDS_CN
        do ib = istartaerband_cn(n),iendaerband_cn(n)
          planckaerband_cn(n) = planckaerband_cn(n) + src1nb(ib)
        end do
      end do
 
!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the non-
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_fr(:,:) = 0.0
      nivl1aer_fr(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_fr(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_fr(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_fr(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_fr(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_FR ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_fr(nband) = nivl + 1
            else
              nivl1aer_fr(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband .eq. 1 .and.   &
              iendsfbands(nivl-1) >= istartaerband_fr(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_fr(1)) then
            nivl1aer_fr(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if (nw >= iendaerband_fr(N_AEROSOL_BANDS_FR) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the 
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_co(:,:) = 0.0
      nivl1aer_co(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_co(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_co(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_co(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_co(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_CO ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_co(nband) = nivl + 1
            else
              nivl1aer_co(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.  &
              iendsfbands(nivl-1) >= istartaerband_co(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_co(1)) then
            nivl1aer_co(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if ( nw >= iendaerband_co(N_AEROSOL_BANDS_CO) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the 
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_cn(:,:) = 0.0
      nivl1aer_cn(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_cn(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_cn(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_cn(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_cn(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_CN ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_cn(nband) = nivl + 1
            else
              nivl1aer_cn(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.  &
              iendsfbands(nivl-1) >= istartaerband_cn(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_cn(1)) then
            nivl1aer_cn(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if ( nw >= iendaerband_cn(N_AEROSOL_BANDS_CN) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the planck-function-weighted band weights for the aerosol
!    parameterization bands onto the non-continuum and continuum ir 
!    aerosol emissivity bands.
!--------------------------------------------------------------------
      sflwwts_fr(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_FR
        do ni=nivl1aer_fr(n),nivl2aer_fr(n)
          sflwwts_fr(n,ni) = planckivlaer_fr(n,ni)/planckaerband(n)
        end do
      end do
      sflwwts_co(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_CO
        do ni=nivl1aer_co(n),nivl2aer_co(n)
          sflwwts_co(n,ni) = planckivlaer_co(n,ni)/     &
                             planckaerband(N_AEROSOL_BANDS_FR+n)
        end do
      end do
      sflwwts_cn(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_CN
        do ni=nivl1aer_cn(n),nivl2aer_cn(n)
          sflwwts_cn(n,ni) = planckivlaer_cn(n,ni)/     &
                             planckaerband_cn(n)
        end do
      end do

!--------------------------------------------------------------------
!    consolidate the continuum and non-continuum weights into an
!    array covering all ir aerosol emissivity bands.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_fr(n,ni)
        end do
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_co(n-N_AEROSOL_BANDS_FR,ni)
        end do
      end do

!-----------------------------------------------------------------
!    allocate and initialize the arrays that will contain the
!    ir aerosol properties for each aerosol optical type
!    over each ir aerosol emissivity band.
!----------------------------------------------------------------
      allocate     &
         (aerextbandlw_MOD   (N_AEROSOL_BANDS, naermodels), &
          aerssalbbandlw_MOD (N_AEROSOL_BANDS, naermodels), &
          aerextbandlw_cn_MOD   (N_AEROSOL_BANDS_CN, naermodels), &
          aerssalbbandlw_cn_MOD (N_AEROSOL_BANDS_CN, naermodels) )
      aerextbandlw_MOD   = 0.0E+00
      aerssalbbandlw_MOD = 0.0E+00
      aerextbandlw_cn_MOD   = 0.0E+00
      aerssalbbandlw_cn_MOD = 0.0E+00

!---------------------------------------------------------------------
!    if longwave aerosol effects are desired, and the following cal-
!    culation has not already been done, calculate the aerosol 
!    properties for each aerosol properties type nw over each aerosol 
!    emissivity band na using the weighted contributions from each
!    aerosol parameterization band ni. mark the calculation as com-
!    pleted.
!
!    the units of extinction coefficient (aeroextivl) are m**2/gm.
!    to make the lw band extinction coefficient (aerextbandlw) have
!    units (m**2/Kg) consistent with the units in FMS models, one
!    must multiply by 1000. this is done below.
!---------------------------------------------------------------------
    if (.not. force_to_repro_quebec) then
      do nw=1,naermodels    
        do na=1,N_AEROSOL_BANDS  
          do ni=1,num_wavenumbers 
            aerextbandlw_MOD(na,nw) = aerextbandlw_MOD(na,nw) + &
                              aeroextivl(ni,nw)*sflwwts(na,ni)*1.0E+03
            aerssalbbandlw_MOD(na,nw) = aerssalbbandlw_MOD(na,nw) + &
                                   aerossalbivl(ni,nw)*sflwwts(na,ni)
          end do
        end do
      end do
      do nw=1,naermodels    
        do na=1,N_AEROSOL_BANDS_CN
          do ni=1,num_wavenumbers 
            aerextbandlw_cn_MOD(na,nw) = aerextbandlw_cn_MOD(na,nw) +&
                          aeroextivl(ni,nw)*sflwwts_cn(na,ni)*1.0E+03
            aerssalbbandlw_cn_MOD(na,nw) =    &
                             aerssalbbandlw_cn_MOD(na,nw) +&
                                aerossalbivl(ni,nw)*sflwwts_cn(na,ni)
          end do
        end do
      end do
    endif

!----------------------------------------------------------------------

end subroutine lw_aerosol_interaction

!######################################################################
! <SUBROUTINE NAME="lw_optical_depth_aerosol">
!  <OVERVIEW>
!   Subroutine to compute aerosol optical depths
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute aerosol optical depths. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_optical_depth_aerosol (Atmos_input, n, Aerosol,    &
!                                  ivol, Optical)
!  </TEMPLATE>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data to model grid point for radiative 
!   properties calculation
!  </IN>
!  <IN NAME="n" TYPE="integer">
!   aerosol optical index
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatological input data
!  </IN>
!  <INOUT NAME="ivol" TYPE="integer">
!   Indexing for radiative properties
!  </INOUT>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Aerosol Optical depth output
!  </INOUT>
! </SUBROUTINE>
!

subroutine lw_optical_depth_aerosol ( relhum, n, Aerosol,  &    
                                      ivol, Aerosolrad_diags, &
                                      aerooptdep )

!------------------------------------------------------------------
! arguments
real, dimension(:,:,:),      intent(in)    :: relhum
integer,                     intent(in)    :: n 
type(aerosol_type),          intent(in)    :: Aerosol
integer, dimension(:,:,:),   intent(in)    :: ivol
type(aerosolrad_diag_type),  intent(inout) :: Aerosolrad_diags
real, intent(out), dimension (size(Aerosol%aerosol,1),  &
                              size(Aerosol%aerosol,2),  &
                              size(Aerosol%aerosol,3))  :: aerooptdep
!------------------------------------------------------------------
! local variables
      real, dimension (size(Aerosol%aerosol,1),  &
                       size(Aerosol%aerosol,2),  &
                       size(Aerosol%aerosol,3), &
                       size(Aerosol%aerosol,4))   :: aerooptdepspec, &
                                                     aerooptdepspec_cn
      integer, dimension (size(Aerosol%aerosol,1),  &
                          size(Aerosol%aerosol,2),  &
                          size(Aerosol%aerosol,3),  &
                          NUM_AERO_INDICES) :: opt_indices

    ! real, dimension (size(Aerosol%aerosol,3)+1) :: bsum

      real      :: asum
      integer   :: nfields, irh
      integer   ::  N_AEROSOL_BANDS
      integer   :: i,j,k
      integer   :: ix, jx, kx
      integer   :: nsc, opt_index, indx
!--------------------------------------------------------------------
      ix = size (Aerosol%aerosol,1)
      jx = size (Aerosol%aerosol,2)
      kx = size (Aerosol%aerosol,3)
      nfields = size (Aerosol%aerosol,4)
!---------------------------------------------------------------------
      aerooptdep(:,:,:) = 0.0

      do k = 1,kx
        do j = 1,jx
          do i = 1,ix
            irh = MIN(100, MAX(0, NINT(100.*relhum(i,j,k))))
            opt_indices(i,j,k,1) = sulfate_index_MOD (irh, ivol(i,j,k) )
            opt_indices(i,j,k,2) = omphilic_index_MOD( irh )
            opt_indices(i,j,k,3) = bcphilic_index_MOD( irh )
            opt_indices(i,j,k,4) = seasalt1_index_MOD( irh )
            opt_indices(i,j,k,5) = seasalt2_index_MOD( irh )
            opt_indices(i,j,k,6) = seasalt3_index_MOD( irh )
            opt_indices(i,j,k,7) = seasalt4_index_MOD( irh )
            opt_indices(i,j,k,8) = seasalt5_index_MOD( irh )
            opt_indices(i,j,k,9) = seasalt_aitken_index_MOD( irh )
            opt_indices(i,j,k,10) = seasalt_fine_index_MOD( irh )
            opt_indices(i,j,k,11) = seasalt_coarse_index_MOD( irh )
            opt_indices(i,j,k,12) = nitrate_index_MOD( irh )
          end do
        end do
      end do

!---------------------------------------------------------------------
!    using relative humidity criterion (where necessary) determine the
!    aerosol category (as an index) appropriate for the aerosol species
!---------------------------------------------------------------------
  do nsc=1, nfields  ! loop on aerosol species

      if (optical_index_MOD(nsc) > 0 ) then
          opt_index = optical_index_MOD(nsc)
          do k = 1, kx
          do j = 1, jx
          do i = 1, ix
              aerooptdepspec(i,j,k,nsc) =    &
                         diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                         (1.0 - aerssalbbandlw_MOD(n,opt_index))*&
                                aerextbandlw_MOD(n,opt_index)
              if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                       diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                       (1.0 - aerssalbbandlw_cn_MOD(n,opt_index))*&
                              aerextbandlw_cn_MOD(n,opt_index)
              end if
          end do
          end do
          end do

      else  ! optical_index_MOD <= 0
          indx = 1 - optical_index_MOD(nsc)
          if (optical_index_MOD(nsc) == BC_FLAG) indx = 1
          if (optical_index_MOD(nsc) == BCPHILIC_FLAG .and. using_im_bcsul) indx = 1

          if (indx < 1 .or. indx > size(opt_indices,4)) then
              call error_mesg ('lw_optical_depth_aerosol', &
                  'Cannot find aerosol optical properties for species = ' // &
                   TRIM( Aerosol%aerosol_names(nsc) ),  FATAL )
          endif

          do k = 1, kx
          do j = 1, jx
          do i = 1, ix
              opt_index = opt_indices(i,j,k,indx)
              aerooptdepspec(i,j,k,nsc) =     &
                        diffac*Aerosol%aerosol(i,j,k,nsc)*&
                        (1.0 - aerssalbbandlw_MOD(n,opt_index))* &
                               aerextbandlw_MOD(n,opt_index)
              if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                          diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                      (1.0 - aerssalbbandlw_cn_MOD(n,opt_index))*&
                             aerextbandlw_cn_MOD(n,opt_index)
              endif
          end do
          end do
          end do

      endif
  end do ! end loop on aerosol species

!---------------------------------------------------------------------
!    save optical path contributions from each layer for band4 and the
!    continuum band. note that if the lw scheme is changed to allow
!    longwave scattering then the %absopdep must be defined approp-
!    riately.
!---------------------------------------------------------------------
  if (n == 1) then
      Aerosolrad_diags%extopdep(:,:,:,:,3) = aerooptdepspec_cn(:,:,:,:)
      Aerosolrad_diags%absopdep(:,:,:,:,3) = aerooptdepspec_cn(:,:,:,:)
  endif

!---------------------------------------------------------------------
!    sum optical depths over all species and obtain column optical depth
!---------------------------------------------------------------------
  do k=1,kx
  do j=1,jx
  do i=1,ix
      asum = 0.0
      do nsc=1,nfields
          asum = asum + aerooptdepspec(i,j,k,nsc)
      end do
      aerooptdep(i,j,k) = asum
  end do
  end do
  end do

!---------------------------------------------------------------------

end subroutine lw_optical_depth_aerosol

!######################################################################

subroutine lw_volcanic_aerosol ( deltaz, n, Aerosolrad_diags, aerooptdep )

!------------------------------------------------------------------
! arguments
real, dimension (:,:,:),     intent(in)    :: deltaz
integer,                     intent(in)    :: n
type(aerosolrad_diag_type),  intent(inout) :: Aerosolrad_diags
real, intent(out), dimension (size(Aerosolrad_diags%lw_ext,1),  &
                              size(Aerosolrad_diags%lw_ext,2),  &
                              size(Aerosolrad_diags%lw_ext,3))  :: aerooptdep
!------------------------------------------------------------------
! local variables
integer :: i, j, k, ix, jx, kx
!---------------------------------------------------------------------
    ix = size (Aerosolrad_diags%lw_ext,1)
    jx = size (Aerosolrad_diags%lw_ext,2)
    kx = size (Aerosolrad_diags%lw_ext,3)

!BW if (Aerosolrad_control%volcanic_lw_aerosols) then
        if (size(Aerosolrad_diags%lw_ext,4) /= 0) then

            do k=2,kx+1
            do j=1,jx
            do i=1,ix
                aerooptdep(i,j,k-1) = Aerosolrad_diags%lw_ext(i,j,k-1,n)*deltaz(i,j,k-1)
            end do
            end do
            end do

            ! diagnostics
            if (n == 5 .or. n == 6) then
              do k=2,kx+1
              do j=1,jx
              do i=1,ix
                Aerosolrad_diags%lw_extopdep_vlcno(i,j,k,n-4) =  &
                                   Aerosolrad_diags%lw_ext(i,j,k-1,n)*deltaz(i,j,k-1)
                Aerosolrad_diags%lw_absopdep_vlcno(i,j,k,n-4) =  &
                                   Aerosolrad_diags%lw_extopdep_vlcno(i,j,k,n-4)
             !! NOT CURRENTLY AVAILABLE IN SEA LW CODE -- lw_ssa not processed
             !! Aerosolrad_diags%lw_absopdep_vlcno(i,j,k,n-4) =  &
             !!              (1.0-Aerosolrad_diags%lw_ssa(i,j,k-1,n))*  &
             !!                   Aerosolrad_diags%lw_ext(i,j,k-1,n)*deltaz(i,j,k-1)
              end do
              end do
              end do
            endif

        endif ! (size)
!BW endif  ! (volcanic_lw_aerosols)

!---------------------------------------------------------------------

end subroutine lw_volcanic_aerosol

!######################################################################
! <SUBROUTINE NAME="sw_aerosol_optical_props">
!  <OVERVIEW>
!   Subroutine that uses the delta-eddington technique in conjunction
!   with a multi-band parameterization for h2o+co2+o2+o3 absorption
!   in the solar spectrum to derive solar fluxes and heating rates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This subroutine calculates optical depth, single scattering albedo,
!    asymmetry parameter of a layer based on gaseous absorbers,
!    clouds, aerosols, and rayleigh scattering. It then uses delta-
!    eddington technique to calculate radiative flux at each layer. 
!    Doubling and adding technique is used to combine the layers
!    and calculate flux at TOA and surface and heating rate. This
!    subroutine allocates a substantial amount of memory and deallocates
!    the allocated memory at the end of the subroutine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call comput(is, ie, js, je, Atmos_input, Surface, Rad_gases, Aerosol, 
!               Astro, &
!               Cldrad_props, Cld_spec)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!    starting subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!    ending subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!    starting subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!    ending subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    Atmos_input_type variable containing the atmospheric
!    input fields on the radiation grid 
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol input data for shortwave radiation calculation
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!    Astronomy_type variable containing the astronomical
!    input fields on the radiation grid  
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!    Radiative_gases_type variable containing the radiative 
!    gas input fields on the radiation grid 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!    The cloud radiative property input fields on the
!    radiation grid
!  </IN>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface data as boundary condition to radiation
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud specification data as initial condition to radiation
!  </IN>
! </SUBROUTINE>
subroutine sw_aerosol_optical_props    &
          (relhum, deltaz, Aerosol, ivol, including_volcanoes,  &
           Aerosolrad_diags, including_aerosols, do_cmip_diagnostics, &
           daylight, aeroextopdep, aerosctopdep, aeroasymfac)
!          aeroextopdep_novolc, aerosctopdep_novolc, aeroasymfac_novolc)

!----------------------------------------------------------------------
!    comput uses the delta-eddington technique in conjunction with a    
!    multiple-band parameterization for h2o+co2+o2+o3 absorption to   
!    derive solar fluxes and heating rates.                             
!    notes: drops are assumed if temp>273.15K, ice crystals otherwise.
!-------------------------------------------------------------------

real, dimension(:,:,:),        intent(in)    :: relhum, deltaz
type(aerosol_type),            intent(in)    :: Aerosol     
integer, dimension(:,:,:),     intent(in)    :: ivol
logical,                       intent(in)    :: including_volcanoes
type(aerosolrad_diag_type),    intent(inout) :: Aerosolrad_diags
logical,                       intent(in)    :: including_aerosols
logical,                       intent(in)    :: do_cmip_diagnostics
logical,dimension(:,:),        intent(in)    :: daylight
real, dimension(:,:,:,:),      intent(out)   :: aeroextopdep, &
                                                aerosctopdep, &
                                                aeroasymfac 
!real, dimension(:,:,:,:), optional, intent(out) :: aeroextopdep_novolc, &
!                                                   aerosctopdep_novolc, &
!                                                   aeroasymfac_novolc 
!

!-------------------------------------------------------------------
!  intent(in) variables:
!
!      relhum         mean relative humidity in model layers
!      deltaz         thickness of model layers
!      Aerosol        aerosol_type structure, contains variables
!                     defining aerosol fields, passed through to
!                     lower level routines
!      ivol           indexing for aerosol radiative properties
!                                                                 
!   intent(inout) variables:
!
!      Aerosolrad_diags  aerosol diagnostics
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables:
 
      real, dimension (size(relhum,3))  :: &
                      arprod, asymm,   arprod2,        &     
                      sum_g_omega_tau, sum_ext,      sum_sct

      integer :: irh
      integer, dimension (size(relhum,3),NUM_AERO_INDICES) :: opt_indices

      real, dimension (size(aerextband_MOD,2))   ::            &
                      aerext,          aerssalb,       aerasymm

      real        :: aerext_i, aerssalb_i, aerasymm_i
      integer     :: j, i, k, nband, nsc, naerosoltypes_used, indx
      integer     :: israd, jsrad, ierad, jerad, ksrad, kerad

      integer :: w340_band_indx, w380_band_indx, w440_band_indx, &
                 w550_band_indx, w670_band_indx, w870_band_indx, &
                 one_micron_indx
      integer :: nbands

!-----------------------------------------------------------------------
!     local variables:
!
!       aeramt
!       sum_g_omega_tau
!       opt_index_v3
!       irh
!    etc.
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------
      israd = 1
      jsrad = 1
      ksrad = 1
      ierad = size(relhum,1)
      jerad = size(relhum,2)
      kerad = size(relhum,3)

      naerosoltypes_used = size(Aerosol%aerosol,4)
      call esfsw_number_of_bands (nbands)

      ! define band indices corresponding to specific wavelengths
      if (do_cmip_diagnostics) then
         call esfsw_bands ( w340_band_indx, w380_band_indx, &
                            w440_band_indx, w550_band_indx, &
                            w670_band_indx, w870_band_indx, &
                            one_micron_indx )
      endif

      do j = JSRAD,JERAD
        do i = ISRAD,IERAD
          if (daylight(i,j) .or. do_cmip_diagnostics) then

            do nband = 1,nbands
              if (including_aerosols) then                           
                aerext(:) = aerextband_MOD(nband,:)
                aerssalb(:) = aerssalbband_MOD(nband,:)
                aerasymm(:) = aerasymmband_MOD(nband,:)

!-------------------------------------------------------------------
!    define the local variables for the band values of aerosol and 
!    cloud single scattering parameters.                               
!    note: the unit for the aerosol extinction is kilometer**(-1).     
!--------------------------------------------------------------------
                do k = KSRAD,KERAD
                  irh = MIN(100, MAX( 0, NINT(100.*relhum(i,j,k))))
                  opt_indices(k,1) = sulfate_index_MOD (irh, ivol(i,j,k) )
                  opt_indices(k,2) = omphilic_index_MOD( irh )
                  opt_indices(k,3) = bcphilic_index_MOD( irh )
                  opt_indices(k,4) = seasalt1_index_MOD( irh )
                  opt_indices(k,5) = seasalt2_index_MOD( irh )
                  opt_indices(k,6) = seasalt3_index_MOD( irh )
                  opt_indices(k,7) = seasalt4_index_MOD( irh )
                  opt_indices(k,8) = seasalt5_index_MOD( irh )
                  opt_indices(k,9) = seasalt_aitken_index_MOD( irh )
                  opt_indices(k,10) = seasalt_fine_index_MOD( irh )
                  opt_indices(k,11) = seasalt_coarse_index_MOD( irh )
                  opt_indices(k,12) = nitrate_index_MOD( irh )
                end do

!---------------------------------------------------------------------
!    calculate scattering properties for all aerosol constituents 
!    combined.
!---------------------------------------------------------------------
                do k = KSRAD,KERAD
                  sum_g_omega_tau(k) = 0.0
                  sum_ext(k) = 0.
                  sum_sct(k) = 0.
                end do
                do nsc = 1,naerosoltypes_used
                  if (optical_index_MOD(nsc) > 0) then
                    aerext_i = aerext(optical_index_MOD(nsc))
                    aerssalb_i = aerssalb(optical_index_MOD(nsc))
                    aerasymm_i = aerasymm(optical_index_MOD(nsc))
                    do k = KSRAD,KERAD
                      arprod(k) = aerext_i*(1.e3*Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb_i*arprod(k)
                      asymm(k)   = aerasymm_i
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) + aerssalb_i*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) +     &
                                      aerasymm_i*(aerssalb_i*arprod(k))
                    end do

                  else ! optical_index_MOD <= 0
                    indx = 1 - optical_index_MOD(nsc)
                    if (optical_index_MOD(nsc) == BC_FLAG) indx = 1
                    if (optical_index_MOD(nsc) == BCPHILIC_FLAG .and. using_im_bcsul) indx = 1

                    if (indx < 1 .or. indx > size(opt_indices,2)) then
                      call error_mesg ('sw_aerosol_optical_props', &
                             'Cannot find aerosol optical properties for species = ' // &
                              TRIM( Aerosol%aerosol_names(nsc) ),  FATAL )
                    end if

                    do k = KSRAD,KERAD
                      arprod(k) = aerext(opt_indices(k,indx)) * (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb(opt_indices(k,indx))*arprod(k)
                      asymm(k)   = aerasymm(opt_indices(k,indx))
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) + aerssalb(opt_indices(k,indx))*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) +  &
                                           aerasymm(opt_indices(k,indx))*  &
                                          (aerssalb(opt_indices(k,indx))*arprod(k))
                    end do

                  endif

                  if (do_cmip_diagnostics) then
                    if (nband == w550_band_indx) then  ! visible light
                      Aerosolrad_diags%extopdep(i,j,:,nsc,1) = arprod(:)
                      Aerosolrad_diags%absopdep(i,j,:,nsc,1) =    &
                                            arprod(:) - arprod2(:)
                      Aerosolrad_diags%asymdep(i,j,:,nsc,1) = asymm(:)
                    endif
                    if (nband == w870_band_indx) then
                      Aerosolrad_diags%extopdep(i,j,:,nsc,6) = arprod(:)
                      Aerosolrad_diags%absopdep(i,j,:,nsc,6) =    &
                                            arprod(:) - arprod2(:)
                      Aerosolrad_diags%asymdep(i,j,:,nsc,6) = asymm(:)
                    endif
                    if (nband == one_micron_indx) then
                      Aerosolrad_diags%extopdep(i,j,:,nsc,2) = arprod(:)
                      Aerosolrad_diags%absopdep(i,j,:,nsc,2) =    &
                                               arprod(:) - arprod2(:)
                      Aerosolrad_diags%asymdep(i,j,:,nsc,2) = asymm(:)
                    endif
                    if (nband == w340_band_indx) then
                      Aerosolrad_diags%extopdep(i,j,:,nsc,7) = arprod(:)
                      Aerosolrad_diags%absopdep(i,j,:,nsc,7) =    &
                                                arprod(:) - arprod2(:)
                      Aerosolrad_diags%asymdep(i,j,:,nsc,7) = asymm(:)
                    endif
                    if (nband == w380_band_indx) then
                      Aerosolrad_diags%extopdep(i,j,:,nsc,8) = arprod(:)
                      Aerosolrad_diags%absopdep(i,j,:,nsc,8) =    &
                                                arprod(:) - arprod2(:)
                      Aerosolrad_diags%asymdep(i,j,:,nsc,8) = asymm(:)
                    endif
                    if (nband == w440_band_indx) then
                      Aerosolrad_diags%extopdep(i,j,:,nsc,9) = arprod(:)
                      Aerosolrad_diags%absopdep(i,j,:,nsc,9) =    &
                                               arprod(:) - arprod2(:)
                      Aerosolrad_diags%asymdep(i,j,:,nsc,9) = asymm(:)
                    endif
                    if (nband == w670_band_indx) then  ! red light
                      Aerosolrad_diags%extopdep(i,j,:,nsc,10) = arprod(:)
                      Aerosolrad_diags%absopdep(i,j,:,nsc,10) =    &
                                                arprod(:) - arprod2(:)
                      Aerosolrad_diags%asymdep(i,j,:,nsc,10) = asymm(:)
                    endif
                  endif
                end do

!----------------------------------------------------------------------
!    if optical depths without volcanic aerosols are requested
!    save them here.
!----------------------------------------------------------------------

               !if (present(aeroextopdep_novolc)) then
               !  do k = KSRAD,KERAD
               !    aeroextopdep_novolc(i,j,k,nband) = sum_ext(k) 
               !    aerosctopdep_novolc(i,j,k,nband) = sum_sct(k) 
               !    aeroasymfac_novolc(i,j,k,nband) = sum_g_omega_tau(k) / (sum_sct(k) + 1.0e-30 )
               !  end do
               !end if

!----------------------------------------------------------------------
!    add the effects of volcanic aerosols, if they are to be included.
!    include generation of diagnostics in the visible (0.55 micron) and
!    nir band (1.0 micron).
!----------------------------------------------------------------------
                if (including_volcanoes) then
                  do k = KSRAD,KERAD
                    sum_ext(k) = sum_ext(k) +    &
                                 Aerosolrad_diags%sw_ext(i,j,k,nband)*  &
                                 deltaz(i,j,k)
                    sum_sct(k) = sum_sct(k) +    &
                                 Aerosolrad_diags%sw_ssa(i,j,k,nband)*  &
                                 Aerosolrad_diags%sw_ext(i,j,k,nband)*  &
                                 deltaz(i,j,k)
                    sum_g_omega_tau(k) =   &
                                 sum_g_omega_tau(k) +&
                                 Aerosolrad_diags%sw_asy(i,j,k,nband)* &
                                 Aerosolrad_diags%sw_ssa(i,j,k,nband)*  &
                                 Aerosolrad_diags%sw_ext(i,j,k,nband)*  &
                                 deltaz(i,j,k)
                    if (do_cmip_diagnostics) then
                      if (nband == w550_band_indx) then   ! visible band
                           Aerosolrad_diags%extopdep_vlcno(i,j,k,1) =   &
                                 Aerosolrad_diags%sw_ext(i,j,k,nband)*  &
                                 deltaz(i,j,k)
                           Aerosolrad_diags%absopdep_vlcno(i,j,k,1) =   &
                            (1.0 - Aerosolrad_diags%sw_ssa(i,j,k,nband))*&
                                Aerosolrad_diags%sw_ext(i,j,k,nband)*  &
                                deltaz(i,j,k)
                      endif
                      if (nband == w870_band_indx) then
                           Aerosolrad_diags%extopdep_vlcno(i,j,k,3) =   &
                                 Aerosolrad_diags%sw_ext(i,j,k,nband)*  &
                                 deltaz(i,j,k)
                           Aerosolrad_diags%absopdep_vlcno(i,j,k,3) =   &
                            (1.0 - Aerosolrad_diags%sw_ssa(i,j,k,nband))*&
                                Aerosolrad_diags%sw_ext(i,j,k,nband)*  &
                                deltaz(i,j,k)
                      endif
                      if (nband == one_micron_indx) then
                           Aerosolrad_diags%extopdep_vlcno(i,j,k,2) =   &
                                 Aerosolrad_diags%sw_ext(i,j,k,nband)*  &
                                 deltaz(i,j,k)
                           Aerosolrad_diags%absopdep_vlcno(i,j,k,2) =   &
                            (1.0 - Aerosolrad_diags%sw_ssa(i,j,k,nband))*&
                                Aerosolrad_diags%sw_ext(i,j,k,nband)*  &
                                deltaz(i,j,k)
                      endif
                    endif
                  end do
                endif   ! (including_volcanoes)
!
!----------------------------------------------------------------------
                do k = KSRAD,KERAD
                  aeroextopdep(i,j,k,nband) = sum_ext(k) 
                  aerosctopdep(i,j,k,nband) = sum_sct(k) 
                  aeroasymfac(i,j,k,nband) = sum_g_omega_tau(k) / (sum_sct(k) + 1.0e-30 )
                end do
              else  ! (if not including_aerosols)
                do k = KSRAD,KERAD
                  aeroextopdep(i,j,k,nband) = 0.0                    
                  aerosctopdep(i,j,k,nband) = 0.0                  
                  aeroasymfac(i,j,k,nband) = 0.0                 
                end do
               !if (present(aeroextopdep_novolc)) then
               !  do k = KSRAD,KERAD
               !    aeroextopdep_novolc(i,j,k,nband) = 0.0
               !    aerosctopdep_novolc(i,j,k,nband) = 0.0
               !    aeroasymfac_novolc(i,j,k,nband) = 0.0
               !  end do
               !end if
              endif ! (including_aerosols)
            end do ! (nband)
          endif  ! (daylight or cmip_diagnostics)

        end do ! (i loop)
      end do   ! (j loop)

!---------------------------------------------------------------------

end subroutine sw_aerosol_optical_props

!#################################################################


               end module aerosolrad_package_mod


