                     module esfsw_driver_mod
!
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
!
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
!
! <OVERVIEW>
!  Code that initializes and calculates shortwave radiative quantities
!  such as flux and heating rate.
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes the necessary shortwave radiative parameters
!  in the initialization subroutine. It then uses delta-eddington approximation
!  and doubling and adding technique to calculate solar flux and
!  heating rate.
! </DESCRIPTION>
!

!    shared modules:

use mpp_mod,              only:  input_nml_file, mpp_chksum
use fms_mod,              only:  open_namelist_file, fms_init, &
                                 mpp_pe, mpp_root_pe, stdlog, &
                                 file_exist, write_version_number, &
                                 check_nml_error, error_mesg, &
                                 FATAL, NOTE, close_file, string, stdout
use constants_mod,        only:  PI, GRAV, radcon_mks, o2mixrat, &
                                 rhoair, pstd_mks, WTMAIR, &
                                 constants_init

!  shared radiation package modules:

use esfsw_parameters_mod, only:  esfsw_parameters_init, &
                                 esfsw_parameters_end

use esfsw_bands_mod,      only:  esfsw_bands_init, &
                                 esfsw_bands_end, &
                                 esfsw_bands, &
                                 esfsw_number_of_bands, &
                                 esfsw_band_segments, &
                                 esfsw_solar_flux, &
                                 esfsw_thickavg

use esfsw_utilities_mod,  only:  esfsw_utilities_init, &
                                 esfsw_utilities_end

use shortwave_types_mod,  only:  sw_output_type

!---------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    esfsw_driver_mod is the internal driver for the esf shortwave
!    package.
!------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public    &
         esfsw_driver_init, swresf,   &
         esfsw_driver_end

private     &

!   called from swresf:
         adding, deledd

!---------------------------------------------------------------------
!  public interfaces/derived-types/assignments
!       inherited from other modules

public  &
  !  from esfsw_bands_mod
        esfsw_number_of_bands, &
        esfsw_bands, &
        esfsw_band_segments, &
        esfsw_solar_flux, &
        esfsw_thickavg

!---------------------------------------------------------------------
!-------- namelist  ---------

logical      ::  do_rayleigh_all_bands = .true. ! rayleigh scattering 
                                                ! calculated in all sw 
                                                ! bands ?
logical      ::  do_herzberg = .false.          ! include the herzberg 
                                                ! effect on the o2 
                                                ! optical depth ?
logical      ::  do_quench = .false.            ! include the quenching
                                                ! effect of non-LTE 
                                                ! processes on the co2 
                                                ! optical depth ?
logical      ::  do_h2o_sw_effects = .true.     ! the shortwave effects
                                                ! of h2o are included ?
logical      ::  do_o3_sw_effects =  .true.     ! the shortwave effects
                                                ! of o3 are included ?
logical      ::  do_co2_sw_effects = .true.     ! the shortwave effects
                                                ! of co2 are included ?                                         
logical      ::  do_o2_sw_effects = .true.      ! the shortwave effects
                                                ! of o2 are included ?
logical      ::  do_ch4_sw_effects = .false.    ! the shortwave effects
                                                ! of ch4 are included ?
logical      ::  do_n2o_sw_effects = .false.    ! the shortwave effects
                                                ! of n2o are included ?
logical      ::  do_sw_continuum =  .false.     ! include the shortwave
                                                ! continuum?
logical      ::  do_rayleigh     = .true.       ! is rayleigh scattering
                                                ! turned on?
logical      ::  reproduce_ulm   = .true.       ! reproduce ulm code
logical      ::  do_four_stream  = .false.      !
                             
logical      ::  remain_rayleigh_bug   = .true.
 
namelist / esfsw_driver_nml /    &
                               do_rayleigh_all_bands, &
                               do_herzberg, do_quench, &
                               do_h2o_sw_effects, do_o3_sw_effects, &
                               do_ch4_sw_effects, do_n2o_sw_effects, &
                               do_co2_sw_effects, do_o2_sw_effects, &
                               do_sw_continuum, do_rayleigh, reproduce_ulm, do_four_stream, remain_rayleigh_bug

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


!---------------------------------------------------------------------
!    variables associated with absorptivity and sw transmission for
!    various gaseous atmospheric components
!
! powph2o      = the scaling factor used in the fit of the h2o          
!                transmission function                                  
!                                                                       
! p0h2o        = the reference pressure (mb) used in the fit of the     
!                h2o transmission function                              
!                                                                       
! c(n)co2(str) = coefficients for the absorptivity expression for co2   
!                for the pressure-scaled and non-scaled, respectively,  
!                portions of the fit (=1.0e-99 if no absorption)        
!                                                                       
! c(n)o2(str)  = coefficients for the absorptivity expression for o2    
!                for the pressure-scaled and non-scaled, respectively,  
!                portions of the fit (=1.0e-99 if no absorption)        
!                                                                       
! ""(schrun)   = coefficients for the absorptivity expression for the   
!                Schuman-Runge o2 band (non-scaled only)                
!                                                                       
! kh2o         =  the psuedo-absorption coefficients in cm2/gm for h2o  
!                                                                       
! ko3          =  the absorption coefficients in cm2/gm for o3    
! kctms         = the psuedo-absorption coefficients in cm2/(molecules*atm)
!                  for h2o self continuum
! kctmf         = the psuedo-absorption coefficients in cm2/(molecules*atm)
!                 for h2o foriegn continuum		     
!                                                                       
! wtfreq       = the weight associated with each exponential term       
! strterm      = logical flag to indicate whether or not a h2o pseudo-  
!                absorption coefficient is assigned a non-scaled        
!                (true) or pressure-scaled (false) gas amount    
!---------------------------------------------------------------------
real, dimension (:), allocatable    :: c1co2, c1co2str, c1o2, c1o2str, &
                                       c2co2, c2co2str, c2o2, c2o2str, &
                                       c3co2, c3co2str, c3o2, c3o2str, &
                                       c4co2, c4co2str, c4o2, c4o2str, &
                                       c1ch4, c1ch4str, c2ch4,         &
                                       c2ch4str, c3ch4, c3ch4str,      &
                                       c4ch4, c4ch4str,                &
                                       c1n2o, c1n2ostr, c2n2o,         &
                                       c2n2ostr, c3n2o, c3n2ostr,      &
                                       c4n2o, c4n2ostr,                &
                                       powph2o, p0h2o
real                                :: c1o2strschrun, c2o2strschrun, &
                                       c3o2strschrun, c4o2strschrun
real, dimension (:), allocatable    :: kh2o, ko3, kctms, kctmf, wtfreq
logical, dimension(:), allocatable  :: strterm
real, dimension (:), allocatable    :: powpco2, p0co2
real, dimension (:,:), allocatable  :: kco2, wtfreqco2
integer, dimension (:), allocatable :: nfreqptsco2
logical, dimension(:,:), allocatable:: strtermco2

!---------------------------------------------------------------------
!    quantities associated with solar spectral parameterization
!                                                                       
! firstrayband = the first band number where the contribution by        
!                rayleigh scattering is included in the solar           
!                calculations                                           
!                                                                       
! nirbands     = the number of bands in the near-infrared (used in      
!                assigning the value of the surface albedo for the      
!                near-infrared, and the visible and ultraviolet         
!                regions, separately)                                   
! nfreqpts     = the number of pseudo-monochromatic frequencies         
! solflxband   = the solar flux in each parameterization band           
! solflxbandref = the solar flux in each parameterization band, used for
!                 defining band-averaged optical parameters. If the
!                 solar constant is time-invariant, it is also the solar
!                 flux in each parameterization band (solflxband).
!---------------------------------------------------------------------
real                                :: refquanray, solflxtotal
integer                             :: firstrayband, nirbands
integer, dimension (:), allocatable :: nfreqpts
real,    dimension(:), allocatable  :: solflxband
real, dimension(:), allocatable     :: wtstr, cosangstr
real, dimension(4)                  :: wtstr_4 =      &
                                       (/0.347854845, 0.652145155,&
                                         0.347854845, 0.652145155/)

integer :: nbands, tot_wvnums, nfrqpts, nh2obands, nstreams
logical :: nstr4 = .false.
integer :: onepsix_band_indx
integer :: visible_band_indx
integer :: nco2bands

!---------------------------------------------------------------------
!    variables associated with rayleigh scattering
!---------------------------------------------------------------------
real, dimension (:), allocatable    :: betaddensitymol

!----------------------------------------------------------------------
!    variables associated with total optical path of species ? - smf
!----------------------------------------------------------------------
real                            :: toto2strmaxschrun
real, dimension(:), allocatable :: totco2max, totco2strmax, &
                                   toto2max, toto2strmax, &
                                   totch4max, totch4strmax, &
                                   totn2omax, totn2ostrmax

!----------------------------------------------------------------------
!    variables associated with the herzberg effect. wtmo2 is the mol-
!    ecular weight of o2. herzberg_fac is a factor used in the last 
!    shortwave band, so that herzberg_fac*wo2 yields a correction for 
!    the o2 optical depth to account for the herzberg o2 heating. this 
!    is done only when do_herzberg is true.
!----------------------------------------------------------------------
real, parameter   :: wtmo2        = 3.19988E+01  
real, parameter   :: wtmco2       = 4.40098E+01
real, parameter   :: herzberg_fac = 9.9488377E-3

!----------------------------------------------------------------------
!    variables associated with the quenching effect. co2_quenchfac is a
!    multiplication factor that reduces the co2 gas optical depth, and 
!    hence, the solar heating in the upper atmosphere, to account for 
!    "quenching" due to non-LTE processes. co2_quenchfac_height is the 
!    reference height for co2_quenchfac [ meters ].
!----------------------------------------------------------------------
real, dimension(30) :: co2_quenchfac
data co2_quenchfac /1.0,.954,.909,.853,.800,.747,.693,.637,.583, .526,&
                    .467,.416,.368,.325,.285,.253,.229,.206,.186,.170,&
                    .163,.156,.151,.144,.138,.132,.127,.124,.068,.037/

real, dimension(30) :: co2_quenchfac_height
data co2_quenchfac_height /67304.,68310.,69303.,70288.,71267.,72245.,&
                           73221.,74195.,75169.,76141.,77112.,78082.,&
                           79051.,80018.,80985.,81950.,82914.,83876.,&
                           84838.,85798.,86757.,87715.,88672.,89627.,&
                           90582.,91535.,92487.,93438.,94387.,106747./



!---------------------------------------------------------------------
!    miscellaneous variables
!---------------------------------------------------------------------
integer, parameter :: NSOLWG = 1
real, dimension(NSOLWG) :: gausswt
logical        :: module_is_initialized = .false.
logical        :: do_esfsw_band_diagnostics = .false.


!---------------------------------------------------------------------
!---------------------------------------------------------------------
 


                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
! <SUBROUTINE NAME="esfsw_driver_init">
!  <OVERVIEW>
!   Subroutine that defines the time-independent quantities associated
!   with the incoming shortwave radiation in the multiple-band solar
!   radiation parameterization.
!  </OVERVIEW>
!  <DESCRIPTION>
!   It first reads in the input namelist and then allocates gas absorption
!   coefficient variables. It then reads in the shortwave input namelist
!   file and assigns the gas absorption coefficients. Rayleigh scattering
!   coefficient is also calculated based on the temperature and pressure
!   structure of the atmosphere.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call esfsw_driver_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine esfsw_driver_init
 
!---------------------------------------------------------------------- 
!    esfsw_driver_init is the constructor for esfsw_driver_mod. it
!    defines the time-independent quantities associated with the 
!    incoming shortwave radiation in the multiple-band solar radiation 
!    parameterization.                                    
!---------------------------------------------------------------------- 

!---------------------------------------------------------------------
!  local variables:

      integer, dimension(:), allocatable  :: endwvnbands
      real,    dimension(:), allocatable  :: solflxbandref
      real,    dimension(:), allocatable  :: solarfluxtoa

      real,    dimension(:), allocatable  :: freqnu
      real,    dimension(:), allocatable  :: ptstr 

      integer, dimension(:), allocatable  :: nwvnsolar
      real   , dimension(:), allocatable  :: solint   

      character(len=64)    :: file_name
      real, dimension(4)   :: ptstr_4 = (/-0.861136312,&
                                          -0.339981044, &
                                           0.861136312,  &
                                           0.339981044 /)
      real      :: ptstr_1 = 0.2
      real      :: temprefray  = 288.15
      real      :: pressrefray = 101325.       ! MKS units
      real      :: densmolref  = 2.54743E+19
      real      :: convfac     = 1.0E+18
      real      :: corrfac, gamma, f1, f2, f3, pinteg, &
                   twopiesq, densmolrefsqt3, wavelength,  &
                   freqsq, ristdm1, ri
      integer   :: iounit, nband, nf, ni, nw, nw1, nw2, nintsolar
      integer   :: unit, io, ierr, logunit
      integer   :: i
      integer   :: n
      real      :: input_flag = 1.0e-99
      
!---------------------------------------------------------------------
!  local variables:
!                                                                       
!      freqnu   
!      ptstr          gaussian points and weights for evaluation of the
!                     diffuse beam.
!      nwvnsolar      the number of wavenumbers in each region where the
!                     solar flux is constant                         
!      solint         the solar flux in watts per meter**2 in each      
!                     wavenumber region where it is constant    
!      endwvnbands    the wavenumber value for the band limits   
!      file_name
!      ptstr_4    
!      ptstr_1     
!      temprefray     reference temperature used in defining rayleigh
!                     optical depth
!      pressrefray    reference pressure used in defining rayleigh
!                     optical depth [ Pa ]
!      densmolref     reference density used in defining rayleigh
!                     optical depth
!      convfac     
!      corrfac
!      gamma
!      f1
!      f2
!      f3
!      pinteg
!      twopiesq
!      densmolrefsqt3
!      wavelength
!      freqsq
!      ristdm1
!      ri
!      iounit
!      nband
!      nf
!      ni
!      nw
!      nw1
!      nw2
!      nintsolar      the number of wavenumber regions where the  
!                     solar flux is constant.   
!      unit
!      io
!      ierr
!      i
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
      call esfsw_utilities_init
      call esfsw_parameters_init (nbands, nfrqpts, &
                                  nh2obands, nstreams, tot_wvnums)

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
       read (input_nml_file, nml=esfsw_driver_nml, iostat=io)
       ierr = check_nml_error(io,'esfsw_driver_nml')
#else   
       if ( file_exist('input.nml')) then
         unit =  open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
         read  (unit, nml=esfsw_driver_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'esfsw_driver_nml')
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
                          write (logunit, nml=esfsw_driver_nml)

!---------------------------------------------------------------------
!    define flag indicating if ICA calculations being done.
!    THIS NEEDS TO BE DEFINED ELSEWHERE
!---------------------------------------------------------------------
!DEL  Cldrad_control%do_ica_calcs = do_ica_calcs
!DEL  Cldrad_control%do_ica_calcs_iz = .true.

      allocate ( freqnu(nbands), &
                 ptstr (nstreams) )

      allocate ( endwvnbands  (0:nbands) )
      allocate ( solflxbandref  (nbands) )
      allocate ( solarfluxtoa   (tot_wvnums) )

!---------------------------------------------------------------------
!    allocate module variables
!---------------------------------------------------------------------
      allocate ( betaddensitymol (nbands) )
      allocate ( c1co2   (nh2obands),  &
                 c1co2str(nh2obands),  &
                 c1o2    (nh2obands),  &
                 c1o2str (nh2obands),  &
                 c2co2   (nh2obands),  &
                 c2co2str(nh2obands),  &
                 c2o2    (nh2obands),  &
                 c2o2str (nh2obands),  &
                 c3co2   (nh2obands),  &
                 c3co2str(nh2obands),  &
                 c3o2    (nh2obands),  &
                 c3o2str (nh2obands),  &
                 c4co2   (nh2obands),  &
                 c4co2str(nh2obands),  &
                 c4o2    (nh2obands),  &
                 c4o2str (nh2obands),  &
                 c1ch4   (nh2obands),  &
                 c1ch4str(nh2obands),  &
                 c2ch4   (nh2obands),  &
                 c2ch4str(nh2obands),  &
                 c3ch4   (nh2obands),  &
                 c3ch4str(nh2obands),  &
                 c4ch4   (nh2obands),  &
                 c4ch4str(nh2obands),  &
                 c1n2o   (nh2obands),  &
                 c1n2ostr(nh2obands),  &
                 c2n2o   (nh2obands),  &
                 c2n2ostr(nh2obands),  &
                 c3n2o   (nh2obands),  &
                 c3n2ostr(nh2obands),  &
                 c4n2o   (nh2obands),  &
                 c4n2ostr(nh2obands),  &
                 powph2o (nh2obands),  &
                 p0h2o   (nh2obands)    )
      allocate ( nfreqpts        (nbands) )
      allocate ( solflxband      (nbands) )
      allocate ( kh2o            (nfrqpts),   & 
                 ko3             (nfrqpts),   &
                 kctms           (nfrqpts),   & 
                 kctmf           (nfrqpts),   & 
                 wtfreq          (nfrqpts),   &
                 strterm         (nfrqpts)   )
      allocate ( wtstr           (nstreams),   & 
                 cosangstr       (nstreams)  )
      allocate ( totco2max    (nh2obands),     &
                 totco2strmax (nh2obands),     &
                 toto2max     (nh2obands),     &
                 toto2strmax  (nh2obands),   &
                 totch4max    (nh2obands),   &
                 totch4strmax (nh2obands),   &
                 totn2omax    (nh2obands),   &
                 totn2ostrmax (nh2obands)    )
      allocate ( nfreqptsco2 (nh2obands) )
      allocate ( wtfreqco2   (nfrqpts,nh2obands) )
      allocate ( kco2        (nfrqpts,nh2obands) )
      allocate ( strtermco2  (nfrqpts,nh2obands) )
      allocate ( powpco2     (nh2obands) )
      allocate ( p0co2       (nh2obands) )

      betaddensitymol = 0.0 ; c1co2    = 0.0 ; c1co2str = 0.0
      c1o2     = 0.0 ; c1o2str  = 0.0 ; c2co2    = 0.0
      c2co2str = 0.0 ; c2o2     = 0.0 ; c2o2str  = 0.0
      c3co2    = 0.0 ; c3co2str = 0.0 ; c3o2     = 0.0
      c3o2str  = 0.0 ; c4co2    = 0.0 ; c4co2str = 0.0
      c4o2     = 0.0 ; c4o2str  = 0.0 ; powph2o  = 0.0
      p0h2o    = 0.0 ; nfreqpts        = 0.0 ; solflxband      = 0.0
      solflxbandref   = 0.0 ; kh2o            = 0.0 ; ko3             = 0.0
      wtfreq          = 0.0 ; strterm     = .FALSE. ; wtstr           = 0.0
      cosangstr       = 0.0 ; totco2max     = 0.0 ; totco2strmax  = 0.0
      toto2max      = 0.0 ; toto2strmax   = 0.0
      kctms    = 0.0 ; kctmf    = 0.0
      c1ch4    = 0.0 ; c1ch4str = 0.0; c2ch4    = 0.0
      c2ch4str = 0.0 ; c3ch4    = 0.0 ; c3ch4str = 0.0
      c4ch4    = 0.0 ; c4ch4str = 0.0
      totch4max     = 0.0 ; totch4strmax  = 0.0
      c1n2o    = 0.0 ; c1n2ostr = 0.0; c2n2o    = 0.0
      c2n2ostr = 0.0 ; c3n2o    = 0.0 ; c3n2ostr = 0.0
      c4n2o    = 0.0 ; c4n2ostr = 0.0
      totn2omax     = 0.0 ; totn2ostrmax  = 0.0
      endwvnbands(0) = 0
      nfreqptsco2 = 0.0 ; wtfreqco2 = 0.0 ; kco2 = 0.0 ; strtermco2 = .FALSE.
      p0co2 = 0.0 ; powpco2 = 0.0

!---------------------------------------------------------------------
!    allocate local variables.
!---------------------------------------------------------------------
      if (nstreams == 4) then
        ptstr(:) = ptstr_4(:)
        wtstr(:) = wtstr_4(:)
        nstr4 = .true.
      else if (nstreams == 1) then
        ptstr(1) = ptstr_1
        wtstr(1) = 1.0
        nstr4 = .false.
      endif

!---------------------------------------------------------------------
!    read input file for band positions, solar fluxes and band
!    strengths.
!---------------------------------------------------------------------

      if (nbands == 18 .and. nfrqpts == 38) then
        if (do_sw_continuum) then
           file_name = 'INPUT/esf_sw_input_data_n38b18ctm'
        else
           file_name = 'INPUT/esf_sw_input_data_n38b18'
        endif
      else if (reproduce_ulm) then
        if (do_sw_continuum) then
           file_name = 'INPUT/esf_sw_input_data_n74b18ctm'
        else
           file_name = 'INPUT/esf_sw_input_data_n74b18'
        endif
      else
         file_name = 'INPUT/esf_sw_input_data'
      endif
      
      call error_mesg ( 'esfsw_driver_mod', &
          'reading solar band data from file '//trim(file_name), NOTE)
      
      iounit = open_namelist_file (file_name)
      read(iounit,101) ( solflxbandref(nband), nband=1,NBANDS )
      read(iounit,102) ( nfreqpts(nband), nband=1,NBANDS )
      read(iounit,103) ( endwvnbands(nband), nband=1,NBANDS )
      read(iounit,103) FIRSTRAYBAND,NIRBANDS
      read(iounit,104) ( powph2o(nband), nband=1,NH2OBANDS )
      read(iounit,104) ( p0h2o(nband), nband=1,NH2OBANDS )
      if ((nbands == 18 .and. nfrqpts == 38) .or. reproduce_ulm) then
        read(iounit,105) ( c1co2(nband), nband=1,NH2OBANDS )
        read(iounit,105) ( c1co2str(nband), nband=1,NH2OBANDS )
        read(iounit,105) ( c2co2(nband), nband=1,NH2OBANDS )
        read(iounit,105) ( c2co2str(nband), nband=1,NH2OBANDS )
        read(iounit,105) ( c3co2(nband), nband=1,NH2OBANDS )
        read(iounit,105) ( c3co2str(nband), nband=1,NH2OBANDS )
      else
        read(iounit,102) NCO2BANDS
        read(iounit,102) ( nfreqptsco2(nband), nband=1,NCO2BANDS )
        read(iounit,104) ( powpco2(nband), nband=1,NCO2BANDS )
        read(iounit,104) ( p0co2(nband), nband=1,NCO2BANDS )
      endif

      read(iounit,105) ( c1o2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1o2str(nband), nband=1,NH2OBANDS ),  &
                         c1o2strschrun
      read(iounit,105) ( c2o2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2o2str(nband), nband=1,NH2OBANDS ), &
                         c2o2strschrun
      read(iounit,105) ( c3o2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3o2str(nband), nband=1,NH2OBANDS ),  &
                         c3o2strschrun
      read(iounit,105) ( c1ch4(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1ch4str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2ch4(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2ch4str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3ch4(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3ch4str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1n2o(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1n2ostr(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2n2o(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2n2ostr(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3n2o(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3n2ostr(nband), nband=1,NH2OBANDS )
      
      if ((nbands == 18 .and. nfrqpts == 38) .or. reproduce_ulm) then
      else 
        do nband = 1,NCO2BANDS
        do nf = 1,nfreqptsco2(nband)
          read(iounit,206) wtfreqco2(nf,nband),kco2(nf,nband),strtermco2(nf,nband)
        end do
        end do
   206  format(1p,2e16.6,16x,l16)
      endif 

      if ((nbands == 18 .and. nfrqpts == 38) .or. reproduce_ulm) then
        if (do_sw_continuum) then   
           do nf = 1,nfrqpts
              read(iounit,1066) wtfreq(nf),kh2o(nf),ko3(nf),kctms(nf),kctmf(nf),strterm(nf)    
           end do
       else 
           do nf = 1,nfrqpts
              read(iounit,106) wtfreq(nf),kh2o(nf),ko3(nf),strterm (nf)
           end do
        endif
      else
        do nf = 1,nfrqpts
           read(iounit,1066) wtfreq(nf),kh2o(nf),ko3(nf),kctms(nf),kctmf(nf),strterm(nf)    
        end do 
      endif

      read(iounit,107) nintsolar

      allocate ( nwvnsolar (nintsolar) )
      allocate ( solint    (nintsolar) )
 
      if (do_o3_sw_effects) then
      else 
        ko3(:) = 0 
      endif 

      if (do_h2o_sw_effects) then
      else
        kh2o(:) = 0
      endif
 
      do ni = 1,nintsolar
        read(iounit,107) nwvnsolar (ni),solint(ni)
      end do
 
      call close_file (iounit)
 
      if (tot_wvnums /= endwvnbands(nbands)) then
        call error_mesg ('esfsw_driver_mod', &
         ' tot_wvnums = '//string(tot_wvnums)//&
         '; endwvnbands(nbands) = '//string(endwvnbands(nbands)),NOTE)
        call error_mesg ( 'esfsw_driver_mod', &
         ' inconsistency between highest solar spectrum wavenumber '//&
          'in esfsw_parameters_mod and in esfsw_sriver input file', &
                                                           FATAL)
      endif

!---------------------------------------------------------------------
!    define the wavenumber one solar fluxes.                       
!----------------------------------------------------------------- --
      do ni = 1,nintsolar
        if ( ni.eq.1 ) then
          nw1 = 1
        else
          nw1 = nw1 + nwvnsolar(ni-1)
        end if
        nw2 = nw1 + nwvnsolar(ni) - 1
        do nw = nw1,nw2
          solarfluxtoa(nw) = solint(ni)
        end do
      end do

!---------------------------------------------------------------------

      call esfsw_bands_init ( endwvnbands, solflxbandref, solarfluxtoa)

!---------------------------------------------------------------------
! define the band index corresponding to visible light (0.55 microns)
! and the band index corresponding to 1.61 microns.
!---------------------------------------------------------------------

      call esfsw_bands ( w550_band_indx=visible_band_indx, &
                         onepsix_micron_indx=onepsix_band_indx )

!---------------------------------------------------------------------
      deallocate  (solint    )
      deallocate  (nwvnsolar )
 
!---------------------------------------------------------------------
!    override the input file value of firstrayband if the nml control
!    variables indicates rayleigh effects are to be considered in
!    all bands.
!--------------------------------------------------------------------
      if (do_rayleigh_all_bands)  firstrayband = 1

!----------------------------------------------------------------------
!    convert some values to mks to match model units
!--------------------------------------------------------------------
      p0h2o = 1.0E-2/p0h2o  ! invert, and convert mb to mks
      kh2o = kh2o *1.0E-01   ! cgs to mks
      ko3  = ko3 *1.0E-01    ! cgs to mks
      
      if ((nbands == 18 .and. nfrqpts == 38) .or. reproduce_ulm) then
     else
         p0co2(1:NCO2BANDS) = 1.0E-2/p0co2(1:NCO2BANDS)
      endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

      do n=1,NH2OBANDS

        if ((nbands == 18 .and. nfrqpts == 38) .or. reproduce_ulm) then
          if (do_co2_sw_effects      .and.  &
              c1co2(n) /= input_flag .and.  &
            c2co2(n) /= input_flag .and.  &
            c3co2(n) /= input_flag ) then
          c4co2(n) = c1co2(n) * c2co2(n) ** c3co2(n)
          c4co2str(n) = c1co2str(n) * c2co2str(n) ** c3co2str(n)
          totco2max(n) = ( (1.0/c1co2(n) ) + c2co2(n) ** c3co2(n) ) ** &
                       (1.0/c3co2(n) ) - c2co2(n)
            if (nbands == 18) then
              if ( n /= 4) then
                totco2strmax(n) = ( (1.0/c1co2str(n) ) + c2co2str(n) ** &
                                 c3co2str(n) ) ** (1.0/c3co2str(n) ) - &
                                c2co2str(n)
              else 
                totco2strmax(n) = HUGE (c4o2strschrun) 
              endif
            else
              totco2strmax(n) = ( (1.0/c1co2str(n) ) + c2co2str(n) ** &
                                 c3co2str(n) ) ** (1.0/c3co2str(n) ) - &
                                c2co2str(n)
            endif
          else !nbands
            c4co2(n) = 0.0                              
            c4co2str(n) = 0.0
            totco2max(n) = 0.0                                            
            totco2strmax(n) = 0.0
          endif
       
        endif

        if (do_o2_sw_effects .and. &
            c1o2(n) /= input_flag .and.   &
            c2o2(n) /= input_flag .and.   &
            c3o2(n) /= input_flag ) then
          c4o2(n) = c1o2(n) * c2o2(n) ** c3o2(n)
          toto2max(n) = ( (1.0/c1o2(n) ) + c2o2(n) ** c3o2(n) ) ** &
                          (1.0/c3o2(n) ) - c2o2(n)   
        else
          c4o2(n) = 0.0                              
          toto2max(n) = -1                                            
        endif


        if (reproduce_ulm) then 

          if (do_o2_sw_effects .and. &
              c1o2str(n) /= input_flag .and.   &
              c2o2str(n) /= input_flag .and.   &
              c3o2str(n) /= input_flag ) then
            c4o2str(n) = c1o2str(n) * c2o2str(n) ** c3o2str(n)
            if (c3o2str(n) > 3.3E-3) then
              toto2strmax(n) = ( (1.0/c1o2str(n) ) + c2o2str(n) ** &
                                  c3o2str(n) ) ** (1.0/c3o2str(n) ) - &
                                  c2o2str(n)
            else
              toto2strmax(n) = HUGE (c4o2strschrun) 
            endif
          else
            c4o2str(n) = 0.0
            toto2strmax(n) = 0.0
          endif

        else ! reproduce_ulm

          if (do_o2_sw_effects .and. &
              c1o2str(n) /= input_flag .and.   &
              c2o2str(n) /= input_flag .and.   &
              c3o2str(n) /= input_flag ) then
            c4o2str(n) = c1o2str(n) * c2o2str(n) ** c3o2str(n)
            if (c3o2str(n) > 5E-2 ) then
              toto2strmax(n) = ( (1.0/c1o2str(n) ) + c2o2str(n) ** &
                                c3o2str(n) ) ** (1.0/c3o2str(n) ) - &
                                c2o2str(n)
            else
              toto2strmax(n) = 1E10!HUGE (c4o2strschrun) 
            endif
          else
            c4o2str(n) = 0.0
            toto2strmax(n) = -1.0
          endif

        endif ! reproduce_ulm



        if (do_ch4_sw_effects) then
          if (c1ch4(n) /= input_flag .and.  &
              c2ch4(n) /= input_flag .and.  &
              c3ch4(n) /= input_flag ) then
            c4ch4(n) = c1ch4(n) * c2ch4(n) ** c3ch4(n)
            c4ch4str(n) = c1ch4str(n) * c2ch4str(n) ** c3ch4str(n)
            totch4max(n) = ( (1.0/c1ch4(n) ) + c2ch4(n) ** c3ch4(n) ) ** &
                             (1.0/c3ch4(n) ) - c2ch4(n)
            totch4strmax(n) = ( (1.0/c1ch4str(n) ) + c2ch4str(n) ** &
                              c3ch4str(n) ) ** (1.0/c3ch4str(n) ) - &
                              c2ch4str(n)
          else
            c4ch4(n) = 0.
            c4ch4str(n) = 0.
            totch4max(n) = 0.
            totch4strmax(n) = 0.     
          endif
        endif ! do_ch4_sw_effects

        if (do_n2o_sw_effects) then
          if (c1n2o(n) /= input_flag .and.  &
              c2n2o(n) /= input_flag .and.  &
              c3n2o(n) /= input_flag ) then
            c4n2o(n) = c1n2o(n) * c2n2o(n) ** c3n2o(n)
            c4n2ostr(n) = c1n2ostr(n) * c2n2ostr(n) ** c3n2ostr(n)
            totn2omax(n) = ( (1.0/c1n2o(n) ) + c2n2o(n) ** c3n2o(n) ) ** &
                             (1.0/c3n2o(n) ) - c2n2o(n)
            totn2ostrmax(n) = ( (1.0/c1n2ostr(n) ) + c2n2ostr(n) ** &
                                 c3n2ostr(n) ) ** (1.0/c3n2ostr(n) ) - &
                                 c2n2ostr(n)
          else
            c4n2o(n) = 0.
            c4n2ostr(n) = 0.
            totn2omax(n) = 0.
            totn2ostrmax(n) = 0.     
          endif
        endif ! do_ch4_sw_effects

      end do ! NH2OBANDS

      c4o2strschrun = c1o2strschrun * c2o2strschrun ** c3o2strschrun
      toto2strmaxschrun = ( (1.0/c1o2strschrun) + c2o2strschrun ** &
                             c3o2strschrun) ** (1.0/c3o2strschrun) - &
                             c2o2strschrun

!----------------------------------------------------------------------
!    define the wavenumbers to evaluate rayleigh optical depth.      
!-------------------------------------------------------------------
      do nband = 1,NBANDS
        freqnu(nband) = 0.5 * ( endwvnbands(nband-1) +   &
                                endwvnbands(nband) )
      end do
 
!---------------------------------------------------------------------
!    define quantities used to determine the rayleigh optical depth. 
!    notes: refquanray is the quantity which multiplies pressure /  
!           temperature to yield the molecular density.                 
!           betaddensitymol is the quantity which multiples the       
!           molecular density to yield the rayleigh scattering      
!           coefficient.                                           
!           2.79E-02 is the depolorization factor.              
!-----------------------------------------------------------------
      refquanray = densmolref * temprefray / pressrefray 
      if ( remain_rayleigh_bug ) then
        corrfac = ( 6.0E+00 + 3.0E+00 * 1.39E-02 )/( 6.0E+00 - 7.0E+00 * &
                    1.39E-02 )
        gamma = 1.39E-02 / ( 2.0E+00 - 1.39E-02 )
      else
        corrfac = ( 6.0E+00 + 3.0E+00 * 2.79E-02 )/( 6.0E+00 - 7.0E+00 * &
                    2.79E-02 )
        gamma = 2.79E-02 / ( 2.0E+00 - 2.79E-02 )
      endif
      f1 = 7.5E-01 / ( gamma * 2.0E+00 + 1.0E+00 )
      f2 = gamma * 3.0E+00 + 1.0E+00 
      f3 = 1.0E+00 - gamma 
      pinteg = 2.0E+00 * PI * ( 2.0E+00 * f1 * f2 * ( 1.0E+00 + f3 / &
               f2 / 3.0E+00 ) )
      twopiesq = 2.0E+00 *  PI ** 2
      densmolrefsqt3 = 3.0E+00 * densmolref ** 2
 
!---------------------------------------------------------------------
      do nband = 1,NBANDS
        wavelength = 1.0E+04 / freqnu(nband)
        freqsq = 1.0 / ( wavelength ) ** 2
        ristdm1 = ( 6.4328E+03 + 2.94981E+06 / ( 1.46E+02 - freqsq ) + &
                    2.554E+04 / ( 4.1E+01 - freqsq ) ) * 1.0E-08
        ri = ristdm1 + 1
        betaddensitymol(nband) = twopiesq*( ri ** 2 - 1.0E+00 ) ** 2 * &
                                 corrfac / ( densmolrefsqt3 *  &
                                 wavelength ** 4 ) * pinteg * convfac 
      end do
 
      gausswt(1) = 1.0

!---------------------------------------------------------------------
!    define the gaussian angles for evaluation of the diffuse beam.  
!--------------------------------------------------------------------
      do i = 1,nstreams
        cosangstr(i) = ( ptstr(i) + 1. ) * 5.0E-01
      end do

!--------------------------------------------------------------------
!    mark the module as initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------
 101  format( 12f10.4 )
 102  format( 32i4 )
 103  format( 20i6 )
 104  format( 12f10.2 )
 105  format( 1p,16e8.1 )
 1066  format( 1p,5e16.6,l16 )
 106  format( 1p,3e16.6,l16 )
 107  format( i5,1p,e14.5 )
 
!---------------------------------------------------------------------

end subroutine esfsw_driver_init 

!#################################################################
! <SUBROUTINE NAME="swresf">
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
!   call swresf(press, pflux, temp, rh2o, deltaz,
!               albedo_vis_dir, albedo_nir_dir,
!               albedo_vis_dif, albedo_nir_dif,
!               qo3, rrvco2, rrvch4, rrvn2o,
!               rrsun, cosz, fracday,
!               camtsw, cldsct, cldext, cldasymm, Sw_output)
!  </TEMPLATE>
!  <IN NAME="press" TYPE="real">
!    layer pressures (mks units) (does not include sfc)
!  </IN>
!  <IN NAME="pflux" TYPE="real">
!    interface pressures (mks units)
!  </IN>
!  <IN NAME="temp" TYPE="real">
!    layer temperatures (K) (does not include sfc)
!  </IN>
!  <IN NAME="rh2o" TYPE="real">
!    h2o mixing ratio
!  </IN>
!  <IN NAME="deltaz" TYPE="real">
!    layer thickness in meters
!  </IN>
!  <IN NAME="albedo_vis_dir" TYPE="real">
!    UV/visible surface albedo for direct radiation
!  </IN>
!  <IN NAME="albedo_nir_dir" TYPE="real">
!    Near-IR surface albedo for direct radiation
!  </IN>
!  <IN NAME="albedo_vis_dif" TYPE="real">
!    UV/visible surface albedo for diffuse radiation
!  </IN>
!  <IN NAME="albedo_nir_dif" TYPE="real">
!    Near-IR surface albedo for diffuse radiation
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!    Astronomy_type variable containing the astronomical
!    input fields on the radiation grid  
!  </IN>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!    The shortwave radiation calculation result
!  </INOUT>
!  <IN NAME="camtsw" TYPE="real">
!   Cloud amount for shortwave clouds. If stochastic clouds is implemented
!   then cloud amount by shortwave band.
!  </IN>
!  <IN NAME="cldext" TYPE="real">
!   Cloud extinction parameter
!  </IN>
!  <IN NAME="cldsct" TYPE="real">
!   Cloud single scattering albedo
!  </IN>
!  <IN NAME="cldasymm" TYPE="real">
!   Cloud asymmetric parameter
!  </IN>
! </SUBROUTINE>

subroutine swresf (press, pflux, temp, rh2o, deltaz, &
                   albedo_vis_dir, albedo_nir_dir, &
                   albedo_vis_dif, albedo_nir_dif, &
                   qo3, rrvco2, rrvch4, rrvn2o, &
                   solflxband, solar_constant, rrsun, cosz, fracday, &
                   camtsw, cldsct, cldext, cldasymm, &
                   aeroasymfac, aerosctopdep, aeroextopdep, &
                   do_totcld_forcing, flag_stoch, Sw_output)

!----------------------------------------------------------------------
!    swresf uses the delta-eddington technique in conjunction with a    
!    multiple-band parameterization for h2o+co2+o2+o3 absorption to   
!    derive solar fluxes and heating rates.                             
!    notes: drops are assumed if temp>273.15K, ice crystals otherwise.
!-------------------------------------------------------------------

real, dimension(:,:,:),        intent(in)    :: press, pflux, temp, rh2o, deltaz
real, dimension(:,:),          intent(in)    :: albedo_vis_dir, &
                                                albedo_nir_dir, &
                                                albedo_vis_dif, &
                                                albedo_nir_dif
real, dimension(:,:,:),        intent(in)    :: qo3
real,                          intent(in)    :: rrvco2, rrvch4, rrvn2o, rrsun
real, dimension(:,:),          intent(in)    :: cosz, fracday
real, dimension(:,:,:,:),      intent(in)    :: camtsw      
real, dimension(:,:,:,:,:),    intent(in)    :: cldsct, cldext, cldasymm
real,dimension(:,:,:,:),       intent(in)    :: aeroasymfac, aerosctopdep, aeroextopdep
real,                          intent(in)    :: solar_constant
real, dimension(nbands),       intent(in)    :: solflxband
logical,                       intent(in)    :: do_totcld_forcing
integer,                       intent(in)    :: flag_stoch
type(sw_output_type),          intent(inout) :: Sw_output   


!-------------------------------------------------------------------
!  intent(in) variables:
!
!      press          layer pressures (mks units) (does not include sfc)
!      pflux          interface pressures (mks units)
!      temp           layer temperatures (K) (does not include sfc)
!      rh2o           h2o mixing ratio
!      deltaz         layer thickness in meters
!      albedo_vis_dir UV/visible surface albedo for direct radiation
!      albedo_nir_dir Near-IR surface albedo for direct radiation
!      albedo_vis_dif UV/visible surface albedo for diffuse radiation
!      albedo_nir_dif Near-IR surface albedo for diffuse radiation
!      qo3            ozone
!      rrvco2,rrvch4,
!          rrvn2o     radiatively active gases (scalars)
!      Astro          astronomy_type structure
!      camtsw         cloud amount for shortwave clouds,
!                     if stochastic clouds is implemented then
!                     cloud amount by shortwave band
!      cldsct, cldext, cldasymm
!                     cloud radiative properties
!      aeroasymfac, aerosctopdep, aeroextopdep
!                     aerosol radiative properties
!                                                                 
!   intent(inout) variables:
!
!      Sw_output         shortwave radiation output data
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables:
 

      logical, dimension (size(temp,1), &
                          size(temp,2), &
                          size(temp,3))  ::  cloud

      logical, dimension (size(temp,1), &
                          size(temp,2))  ::   &
                              cloud_in_column,   daylight

      real, dimension (size(temp,1), &
                       size(temp,2), &
                       size(temp,3), &
                        nbands)  ::  rayopdep

      real, dimension (size(temp,1), &
                       size(temp,2), &
                       size(temp,3), &
                       nfrqpts,NSOLWG)  ::  gasopdep

      real, dimension (size(temp,1), &
                       size(temp,2), &
                       size(temp,3))  :: &
            cloudfrac,      cldfrac_band,    cldfrac,         &
                    cloud_deltaz,                     &
            cloudasymfac,   cloudextopdep,                     &
            cloudsctopdep,                   deltap,           &
            densitymol,     extopdepclr,                       &
            extopdepovc,    fclr,            fovc,             &
            gocpdp,         gclr,            gstrclr,          &
            gstrovc,        govc,            omegastrclr,      &
            omegastrovc,    rlayerdif,       rlayerdir,        &
            rlayerdifclr,   rlayerdifovc,    &
            rlayerdirclr,   rlayerdirovc,                     &
            sctopdepclr,    sctopdepovc,      &
            ssalbclr,        ssalbovc,                         &
            taustrclr,      taustrovc,        &
            tlayerde,       tlayerdeclr,      &
            tlayerdeovc,    tlayerdif,                        &
            tlayerdifclr,   tlayerdifovc,     &
            tlayerdir,      tlayerdirclr,     &
            tlayerdirovc, pmstr1clr, pmstr2clr, pmstr3clr,  &
            pmstr1ovc, pmstr2ovc, pmstr3ovc
           
      real, dimension (size(pflux,1), &
                       size(pflux,2), &
                       size(pflux,3))  :: &
              reflectance,  transmittance,  tr_dir            

      real, dimension (size(pflux,1), &
                       size(pflux,2), &
                       size(pflux,3)) :: &
            dfswbandclr,     fswbandclr,    ufswbandclr, &
            dfswband,        fswband,       ufswband,  &
            sumtrclr,        sumreclr,      sumtr_dir_clr,  & 
            sumtr,           sumre,         sumtr_dir

      real, dimension (size(temp,1), &
                       size(temp,2)) :: sumtr_dir_up

      real, dimension (size(pflux,1), &
                       size(pflux,2), &
                       size(pflux,3))  :: &
            pflux_cgs,            pflux_mks, &
            reflectanceclr,  transmittanceclr, tr_dirclr

      real, dimension (size(temp,1), &
                       size(temp,2), &
                       size(temp,3)) :: hswbandclr, hswband

      real, dimension (size(temp,1), &
                       size(temp,2))    ::  &
            sfcalb_dir,   sfcalb_dif,  wtfac_p,    &
            fracday_p,    solarflux_p

      real, dimension (size(press,1),    &
                       size(press,2), NSOLWG) :: cosangsolar_p

      real, dimension (size(press,1), &
                       size(press,2), &
                       size(press,3)) :: press_mks

      integer :: j, i, k, ng, np, nband, nf, ns

      integer :: nprofile, nprofiles
      real    :: profiles_inverse

      integer :: ix, jx, kx, israd, jsrad, ierad, jerad, ksrad, kerad
      real    :: ssolar  
      real    :: solflxtotal_local

      logical :: do_stochastic_clouds, do_ica_calcs


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

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('esfsw_driver_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    set flags for stochastic clouds
!---------------------------------------------------------------------
     do_stochastic_clouds = flag_stoch .gt. 0
     do_ica_calcs         = flag_stoch .gt. 1
     if (do_stochastic_clouds) then
        if (size(camtsw,4) .ne. nbands) &
            call error_mesg ('esfsw_driver_mod',   &
               'incorrect size for cloud prop arrays when '//&
               'stochastic clouds activated', FATAL )
     endif
     if (do_ica_calcs) then
        if (size(cldsct,5) .ne. nbands) &
            call error_mesg ('esfsw_driver_mod',   &
               'incorrect size for cloud prop arrays when '//&
               'ica calcs activated', FATAL )
     endif

!---------------------------------------------------------------------
!    define the solar_constant appropriate at Rad_time when solar
!    input is varying with time.
!---------------------------------------------------------------------
        solflxtotal_local = 0.0
        do nband = 1,NBANDS
        solflxtotal_local = solflxtotal_local + solflxband(nband)
        end do

!---------------------------------------------------------------------
!  convert to cgs and then back to mks for consistency with previous 
!---------------------------------------------------------------------
      press_mks(:,:,:) = 0.1*(10.0*press(:,:,:))
      pflux_cgs(:,:,:) =     (10.0*pflux(:,:,:))
      cloud_deltaz = deltaz  ! temporary mod - BW
                             ! until cloud props are updated

!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------
      ix = size(temp,1)
      jx = size(temp,2)
      kx = size(temp,3)
      israd = 1
      jsrad = 1
      ksrad = 1
      ierad = ix
      jerad = jx
      kerad = kx

!----------------------------------------------------------------------c
!    define a flag indicating columns in which there is sunshine during
!    this radiation time step. define a flag indicating points with both
!    sunlight and cloud.      
!----------------------------------------------------------------------c
      do j = JSRAD,JERAD
        do i = ISRAD,IERAD
          if ( fracday(i,j) /= 0.0 ) then
            daylight(i,j) = .true.                 
          else
            daylight(i,j) = .false.                
          endif     
          cloud_in_column(i,j) = .false.
        end do
      end do
        
!---------------------------------------------------------------------

      ssolar = solar_constant*rrsun
 
!----------------------------------------------------------------------
!    define a flag indicating points with both sunlight and cloud. set
!    the flag indicating that cloud is present in columns with such
!    points.
!----------------------------------------------------------------------
      if (.not. do_stochastic_clouds) then
        cldfrac = camtsw(:,:,:,1)
        do j = JSRAD,JERAD
          do i = ISRAD,IERAD
            if (daylight(i,j)) then
              do k=KSRAD,KERAD
                if (cldfrac(i,j,k) > 0.0)  then
                  cloud_in_column(i,j) = .true.
                  cloud(i,j,k) = .true.
                  cloudfrac(i,j,k) = cldfrac(i,j,k)
                else
                  cloud(i,j,k) = .false.
                  cloudfrac(i,j,k) = 0.0
                endif
              end do
            else
              do k=KSRAD,KERAD
                cloud(i,j,k) = .false.
                cloudfrac(i,j,k) = 0.0
              end do
            endif
          end do
        end do
      endif

!----------------------------------------------------------------------c
!    define pressure related quantities, pressure is in mks units. 
!----------------------------------------------------------------------c
      pflux_mks = pflux_cgs*1.0E-1

      do k = KSRAD+1,KERAD+1
        deltap(:,:,k-1) = pflux_mks(:,:,k) - pflux_mks(:,:,k-1)
        gocpdp(:,:,k-1) = radcon_mks/deltap(:,:,k-1)
      end do
 
      call compute_gas_props (press_mks, pflux_mks, temp, rh2o, deltaz, &
                              qo3, rrvco2, rrvch4, rrvn2o, &
                              cosz, daylight, gasopdep)
  
!---------------------------------------------------------------------
!    define the molecular density for use in calculating the           
!    rayleigh optical depth (deltaz is in meters).                     
!--------------------------------------------------------------------
      do k = KSRAD,KERAD
        densitymol(:,:,k) = refquanray * press_mks(:,:,k) / temp(:,:,k)
      end do
 

! assumption is that there is 1 cloud profile for each sw band
      if (do_ica_calcs) then
        nprofiles = nbands
        profiles_inverse = 1.0/nprofiles
      else
        nprofiles = 1
        profiles_inverse = 1.0
      endif

!--------------------------------------------------------------------
!    define the rayleigh optical depths.                                
!---------------------------------------------------------------------
      do nband = 1, NBANDS
        if (do_rayleigh) then
        rayopdep(:,:,:,nband) = betaddensitymol(nband)*  &
                                        densitymol(:,:,:)*deltaz(:,:,:)
        else 
          rayopdep(:,:,:,nband) = 1.e-6*betaddensitymol(nband)*  &
                                          densitymol(:,:,:)*deltaz(:,:,:)
        endif
      end do   ! (nband loop)

!--------------------------------------------------------------------
 
      do nprofile=1, nprofiles
        if (do_ica_calcs) then
          cldfrac_band(:,:,:) = camtsw(:,:,:,nprofile)
        endif

!----------------------------------------------------------------------c
!    np is a counter for the pseudo-monochromatic frequency point 
!    number.
!----------------------------------------------------------------------c
        np = 0
 
        if (do_totcld_forcing) then
          dfswbandclr(:,:,:) = 0.0
          fswbandclr(:,:,:) = 0.0
          hswbandclr(:,:,:) = 0.0
          ufswbandclr(:,:,:) = 0.0
        endif
        reflectanceclr = 0.0
        transmittanceclr = 0.0

!----------------------------------------------------------------------c
!    begin band loop                                                   
!----------------------------------------------------------------------c
        do nband = 1,NBANDS
 
          sumtr(:,:,:) = 0.0
          sumtr_dir(:,:,:) = 0.0
          sumtr_dir_up(:,:) = 0.0
          sumre(:,:,:) = 0.0
          if (do_totcld_forcing) then
            sumtrclr(:,:,:) = 0.0
            sumreclr(:,:,:) = 0.0
            sumtr_dir_clr(:,:,:) = 0.0
          endif
          if (do_stochastic_clouds .and. .not. do_ica_calcs) then
            cldfrac_band(:,:,:) = camtsw(:,:,:,nband)
          endif
 
!----------------------------------------------------------------------
!    if stochastic clouds are activated (cloud fields differ with sw
!    parameterization band), define a flag indicating points with both
!    sunlight and cloud for the current sw parameterization band. set
!    the flag indicating that cloud is present in columns with such
!    points.
!----------------------------------------------------------------------
          if (do_stochastic_clouds) then
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                cloud_in_column(i,j) = .false.
                if (daylight(i,j)) then
                  do k = KSRAD,KERAD
                    if (cldfrac_band(i,j,k) > 0.0)  then
                      cloud_in_column(i,j) = .true.
                      cloud(i,j,k) = .true.
                      cloudfrac(i,j,k) = cldfrac_band(i,j,k)
                    else
                      cloud(i,j,k) = .false.
                      cloudfrac(i,j,k) = 0.0
                    endif
                  end do
                else
                  do k = KSRAD,KERAD
                    cloud(i,j,k) = .false.
                    cloudfrac(i,j,k) = 0.0
                  end do
                endif
              end do
            end do
          endif


!---------------------------------------------------------------------
!    obtain cloud properties from the input arguments
!--------------------------------------------------------------------

          do k = KSRAD,KERAD
            do j=JSRAD,JERAD
              do i=ISRAD, IERAD
                if (cloud(i,j,k) ) then
                  cloudextopdep(i,j,k) = 1.0E-03*cldext(i,j,k,nband,nprofile) * &
                                         cloud_deltaz(i,j,k)
                  cloudsctopdep(i,j,k) = 1.0E-03*cldsct(i,j,k,nband,nprofile) * &
                                         cloud_deltaz(i,j,k)
!DEL              cloudextopdep(i,j,k) = cldext(i,j,k,nband,nprofile)
!DEL              cloudsctopdep(i,j,k) = cldsct(i,j,k,nband,nprofile)
                  cloudasymfac(i,j,k) = cldasymm(i,j,k,nband,nprofile)
                endif
              end do
            end do
          end do

!-----------------------------------------------------------------
!    define clear sky arrays
!-----------------------------------------------------------------
          if (nband >= firstrayband ) then
            do k=ksrad,kerad
              do j=jsrad,jerad
                do i=israd,ierad
                  if (daylight(i,j) ) then
                    sctopdepclr(i,j,k) = rayopdep(i,j,k,nband) +   &
                                            aerosctopdep(i,j,k,nband)
                    gclr(i,j,k) = aeroasymfac(i,j,k,nband)*  &
                                  aerosctopdep(i,j,k,nband)/&
                                                     sctopdepclr(i,j,k)
                    if (do_four_stream) then
                      fclr(i,j,k) = aeroasymfac(i,j,k,nband)**3*gclr(i,j,k)
                    else
                      fclr(i,j,k) = aeroasymfac(i,j,k,nband)*gclr(i,j,k)
                    end if 
                    gstrclr(i,j,k) = ( gclr(i,j,k)  - fclr(i,j,k) )/  &
                                     ( 1.0 - fclr(i,j,k) )
                  endif
                end do
              end do
            end do
          endif

!-----------------------------------------------------------------
!    define cloudy sky arrays
!-----------------------------------------------------------------
          do k=KSRAD,KERAD
            do j=JSRAD,JERAD
              do i=ISRAD,IERAD
                if (cloud(i,j,k)) then
                  sctopdepovc(i,j,k) = rayopdep(i,j,k,nband) +    &
                                       aerosctopdep(i,j,k,nband) + &
                                       cloudsctopdep(i,j,k) 
                  govc(i,j,k) = ( ( cloudasymfac(i,j,k) *   &
                                    cloudsctopdep(i,j,k) ) +  &
                                  ( aeroasymfac(i,j,k,nband) *   &
                                    aerosctopdep(i,j,k,nband)))/   &
                                    sctopdepovc(i,j,k)
                  if (do_four_stream) then
                    fovc(i,j,k) = ( ( cloudasymfac(i,j,k) ** 4 *  &
                                      cloudsctopdep(i,j,k) ) + &
                                    ( aeroasymfac(i,j,k,nband) ** 4 *  &
                                      aerosctopdep(i,j,k,nband) ))/  &
                                      sctopdepovc(i,j,k)
                  else
                    fovc(i,j,k) = ( ( cloudasymfac(i,j,k) ** 2 *  &
                                      cloudsctopdep(i,j,k) ) + &
                                    ( aeroasymfac(i,j,k,nband) ** 2 *  &
                                      aerosctopdep(i,j,k,nband) ))/  &
                                      sctopdepovc(i,j,k)
                  end if 
                  gstrovc(i,j,k) = ( govc(i,j,k)  - fovc(i,j,k))/  &
                                   ( 1.0 - fovc(i,j,k) )
                endif
              end do
            end do
          end do

!---------------------------------------------------------------------
!    begin frequency points in the band loop.                          
!--------------------------------------------------------------------
          do nf = 1,nfreqpts(nband)
            np = np + 1
 
!---------------------------------------------------------------------
!    begin gaussian angle loop (ng > 1 only when lswg = true).        
!--------------------------------------------------------------------
            do ng = 1,NSOLWG
 
!---------------------------------------------------------------------
!    clear sky mode                                                    
!    note: in this mode, the delta-eddington method is performed for all
!    spatial columns experiencing sunlight.         
!--------------------------------------------------------------------
              if (nband >= firstrayband )  then
                do k=ksrad,kerad
                  do j=jsrad,jerad
                    do i=israd,ierad
                      if (daylight(i,j) ) then
                        extopdepclr(i,j,k) = gasopdep(i,j,k,np,ng) +  &
                                             rayopdep(i,j,k,nband) +   &
                                             aeroextopdep(i,j,k,nband)
                        ssalbclr(i,j,k) = sctopdepclr(i,j,k)/    &
                                          extopdepclr(i,j,k)
                        taustrclr(i,j,k) = extopdepclr(i,j,k)*( 1.0 -  &
                                           ssalbclr(i,j,k)*fclr(i,j,k) )
                        omegastrclr(i,j,k) =      &
                               ssalbclr(i,j,k)*((1.0 - fclr(i,j,k))/  &
                               (1.0 -  ssalbclr(i,j,k)*fclr(i,j,k)))

                        !stunew for 4 stream code 
                        pmstr1clr(i,j,k) =        &
                             3.0 *( gclr(i,j,k)-fclr(i,j,k) ) / ( 1.0-fclr(i,j,k) )
                        pmstr2clr(i,j,k) =        &
                             5.0 *( gclr(i,j,k)**2-fclr(i,j,k) ) / ( 1.0-fclr(i,j,k) )
                        pmstr3clr(i,j,k) =        &
                             7.0 *( gclr(i,j,k)**3-fclr(i,j,k) ) / ( 1.0-fclr(i,j,k) )       
                        !stunew for 4 stream code 
                      endif
                    end do
                  end do
                end do
              endif

!--------------------------------------------------------------------
!    calculate the scaled single-scattering quantities for use in the   
!    delta-eddington routine.                                         
!--------------------------------------------------------------------
              do k=KSRAD,KERAD
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    if (cloud(i,j,k) ) then
                      extopdepovc(i,j,k) = gasopdep(i,j,k,np,ng) +    &
                                           rayopdep(i,j,k,nband) +  &
                                           aeroextopdep(i,j,k,nband) + &
                                           cloudextopdep(i,j,k)
                      ssalbovc(i,j,k) = sctopdepovc(i,j,k) /    &
                                        extopdepovc(i,j,k)
                      taustrovc(i,j,k) = extopdepovc(i,j,k)*( 1.0 - &
                                         ssalbovc(i,j,k)*fovc(i,j,k) )
                      omegastrovc(i,j,k) = ssalbovc(i,j,k)*( ( 1.0 - &
                                           fovc(i,j,k) )/( 1.0 -   &
                                           ssalbovc(i,j,k) *   &
                                           fovc(i,j,k) ) )
                      !stunew for 4 stream code                  
                      pmstr1ovc(i,j,k) =        &
                           3.0 *( govc(i,j,k)-fovc(i,j,k) ) / ( 1.0-fovc(i,j,k) )
                      pmstr2ovc(i,j,k) =        &
                           5.0 *( govc(i,j,k)**2-fovc(i,j,k) ) / ( 1.0-fovc(i,j,k) )
                      pmstr3ovc(i,j,k) =        &
                           7.0 *( govc(i,j,k)**3-fovc(i,j,k) ) / ( 1.0-fovc(i,j,k) )
                      !stunew for 4 stream code  
                    endif
                  end do
                end do
              end do

!---------------------------------------------------------------------
!    do calculations for all desired zenith angles.
!---------------------------------------------------------------------
              cosangsolar_p(:,:,ng) = cosz(:,:)   
              where (cosangsolar_p(:,:,:) == 0.0)   &
                                          cosangsolar_p(:,:,:) = 1.0

!---------------------------------------------------------------------
!    clear sky mode                                                    
!    note: in this mode, the delta-eddington method is performed for all
!    spatial columns experiencing sunlight.         
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    calculate the scaled single-scattering quantities for use in the   
!    delta-eddington routine.                                       
!---------------------------------------------------------------------
              if (nband >= firstrayband )  then

!---------------------------------------------------------------------
!    do diffuse calculation only for first zenith angle -- it is 
!    independent of zenith angle.
!---------------------------------------------------------------------
                if (do_four_stream) then
                  call deledd_4stream     &
                      (ix, jx, kx, taustrclr, omegastrclr, &
                       pmstr1clr, pmstr2clr, pmstr3clr, &
                       cosangsolar_p(:,:,ng), ng, daylight, &
                       rlayerdirclr, tlayerdirclr, tlayerdeclr, &
                       rlayerdif=rlayerdifclr, tlayerdif=tlayerdifclr)
                else
                  call deledd     &
                      (ix, jx, kx, taustrclr, omegastrclr, &
                       gstrclr, cosangsolar_p(:,:,ng), ng, daylight, &
                       rlayerdirclr, tlayerdirclr, tlayerdeclr, &
                       rlayerdif=rlayerdifclr, tlayerdif=tlayerdifclr)
                endif

!---------------------------------------------------------------------
!    the following needs to be done at daylight points only -- currently
!    this code is not active, since ng == 1.
!---------------------------------------------------------------------
                if (ng /= 1) then
                  tlayerdifclr = 0.0       
                  if (NSTREAMS == 1) then
                    tlayerdifclr(:,:,:) =   &
                           exp(-gasopdep(:,:,:,np,ng)/cosangstr(1))
                  else
                    do ns = 1,NSTREAMS
                      tlayerdifclr(:,:,:) =    &
                           tlayerdifclr(:,:,:) +  & 
                                 exp( -gasopdep(:,:,:,np,ng)/&
                                       cosangstr(ns) )*wtstr(ns)* &
                                       cosangstr(ns)
                    end do
                  endif
                  rlayerdifclr = 0.0
                endif
  
!---------------------------------------------------------------------
!    initialize the layer reflection and transmission arrays for the   
!    non-rayleigh-scattering case.                            
!-------------------------------------------------------------------
              else
!---------------------------------------------------------------------
!    the following needs to be done at daylight points only -- currently
!    this code is not active, since ng == 1, and all bands see rayleigh
!    scattering.
!---------------------------------------------------------------------
                tlayerdifclr = 0.0       
                if (NSTREAMS == 1) then
                  tlayerdifclr(:,:,:) =     &
                          exp( -gasopdep(:,:,:,np,ng)/cosangstr(1))
                else
                  do ns = 1,NSTREAMS
                    tlayerdifclr(:,:,:) =   &
                       tlayerdifclr(:,:,:) +    &
                            exp(-gasopdep(:,:,:,np,ng)/&
                              cosangstr(ns))*wtstr(ns)*cosangstr(ns)
                  end do
                endif
                rlayerdirclr(:,:,:) = 0.0
                do k=KSRAD,KERAD
                  tlayerdirclr(:,:,k) =      &
                                exp( -gasopdep(:,:,k,np,ng) /   &
                                               cosangsolar_p(:,:,ng) )
                end do  
                tlayerdeclr(:,:,:) = tlayerdirclr(:,:,:)
                rlayerdifclr = 0.0
              endif

!---------------------------------------------------------------------
!    overcast sky mode                                                  
!    note: in this mode, the delta-eddington method is performed only 
!    for spatial columns containing a cloud and experiencing sunlight. 
!---------------------------------------------------------------------


!----------------------------------------------------------------------
!    calculate the reflection and transmission in the scattering layers 
!    using the delta-eddington method.                                  
!-------------------------------------------------------------------
              if (do_four_stream) then
                call deledd_4stream      &
                     (ix, jx, kx, taustrovc, omegastrovc, &
                      pmstr1ovc, pmstr2ovc, pmstr3ovc, &
                      cosangsolar_p(:,:,ng), ng, daylight, &
                      rlayerdirovc, tlayerdirovc, tlayerdeovc, &
                      rlayerdif=rlayerdifovc, tlayerdif=tlayerdifovc,&
                      cloud=cloud)
              else
                call deledd     &
                     (ix, jx, kx, taustrovc, omegastrovc, gstrovc, &
                      cosangsolar_p(:,:,ng), ng, daylight, &
                      rlayerdirovc, tlayerdirovc, tlayerdeovc, &
                      rlayerdif=rlayerdifovc, tlayerdif=tlayerdifovc,&
                      cloud=cloud)
              end if
              if (ng /= 1) then
                tlayerdifovc(:,:,:) = tlayerdifclr(:,:,:)
                rlayerdifovc(:,:,:) = rlayerdifclr(:,:,:)
              endif
 
!-------------------------------------------------------------------- 
!    weight the reflection and transmission arrays for clear and        
!    overcast sky conditions by the cloud fraction, to calculate the    
!    resultant values.                                                  
!---------------------------------------------------------------------- 
              do k=KSRAD,KERAD
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    if ( cloud(i,j,k) ) then
                      rlayerdir(i,j,k) = cloudfrac(i,j,k)*   &
                                         rlayerdirovc(i,j,k) +  &
                                         (1.0 - cloudfrac(i,j,k) )*  &
                                         rlayerdirclr(i,j,k)
                      rlayerdif(i,j,k) = cloudfrac(i,j,k) *  &
                                         rlayerdifovc(i,j,k) +  &
                                         ( 1.0 - cloudfrac(i,j,k) )* &
                                         rlayerdifclr(i,j,k)
                      tlayerdir(i,j,k) = cloudfrac(i,j,k) *   &
                                         tlayerdirovc(i,j,k) +  &
                                         ( 1.0 - cloudfrac(i,j,k) )* &
                                         tlayerdirclr(i,j,k)
                      tlayerdif(i,j,k) = cloudfrac(i,j,k) *   &
                                         tlayerdifovc(i,j,k) +  &
                                         ( 1.0 - cloudfrac(i,j,k) )* &
                                         tlayerdifclr(i,j,k)
                      tlayerde(i,j,k) =  cloudfrac(i,j,k) *   &
                                         tlayerdeovc(i,j,k) +  &
                                         (1.0 - cloudfrac(i,j,k) )* &
                                         tlayerdeclr(i,j,k)
                    else if (daylight(i,j)) then
                      rlayerdir(i,j,k) = rlayerdirclr(i,j,k)
                      tlayerdir(i,j,k) = tlayerdirclr(i,j,k)
                      rlayerdif(i,j,k) = rlayerdifclr(i,j,k)
                      tlayerdif(i,j,k) = tlayerdifclr(i,j,k)
                      tlayerde (i,j,k) = tlayerdeclr (i,j,k)
                    endif
                  end do
                end do
              end do
 
!---------------------------------------------------------------------
!    define the surface albedo (infrared value for infrared bands,      
!    visible value for the remaining bands).                            
!----------------------------------------------------------------------c
              if (nband <= NIRBANDS ) then
                sfcalb_dir(:,:) = albedo_nir_dir(:,:)
                sfcalb_dif(:,:) = albedo_nir_dif(:,:)
              else
                sfcalb_dir(:,:) = albedo_vis_dir(:,:)
                sfcalb_dif(:,:) = albedo_vis_dif(:,:)
              end if
 
!-------------------------------------------------------------------- 
!    calculate the reflection and transmission at flux levels from the  
!    direct and diffuse values of reflection and transmission in the  
!    corresponding layers using the adding method.                      
!---------------------------------------------------------------------
              call adding         &
                  (ix, jx, kx, rlayerdir, tlayerdir, rlayerdif, &
                   tlayerdif, tlayerde, sfcalb_dir, sfcalb_dif,    &
                   daylight, reflectance, transmittance, tr_dir)    

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
              if (do_totcld_forcing) then
                call adding       &
                    (ix, jx,  kx, rlayerdirclr, tlayerdirclr,   &
                     rlayerdifclr, tlayerdifclr, tlayerdeclr,   &
                     sfcalb_dir,  sfcalb_dif, cloud_in_column,  &
                     reflectanceclr, transmittanceclr, tr_dirclr)
              endif

!---------------------------------------------------------------------- 
!    weight and sum the reflectance and transmittance to calculate the 
!    band values.                                                     
!-------------------------------------------------------------------
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  wtfac_p(i,j) = wtfreq(np)*gausswt(ng)*   &
                                                cosangsolar_p(i,j,ng)
                end do
              end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
              do k = KSRAD,KERAD+1
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    if (daylight(i,j) ) then
                      sumtr(i,j,k) = sumtr(i,j,k) +    &
                                    transmittance(i,j,k)*wtfac_p(i,j)
                      sumtr_dir(i,j,k) = sumtr_dir(i,j,k) +  &
                                      tr_dir(i,j,k)*wtfac_p(i,j)
                      sumre(i,j,k) = sumre(i,j,k) +     &
                                     reflectance(i,j,k)*wtfac_p(i,j)
                    endif
                  end do
                end do
              end do

              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  if (daylight(i,j) ) then
                    sumtr_dir_up(i,j) = sumtr_dir_up(i,j) + &
                      tr_dir(i,j,KERAD+1)*sfcalb_dir(i,j)*wtfac_p(i,j)
                  endif
                end do
              end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
              if (do_totcld_forcing) then
                do k = KSRAD,KERAD+1
                  do j=JSRAD,JERAD
                    do i=ISRAD,IERAD
                      if (cloud_in_column(i,j)) then
                        sumtrclr(i,j,k) =    &
                              sumtrclr(i,j,k) +   &
                              transmittanceclr(i,j,k)* wtfac_p(i,j) 
                        sumtr_dir_clr(i,j,k) = &
                              sumtr_dir_clr(i,j,k) +  &
                              tr_dirclr(i,j,k)*wtfac_p(i,j)
                        sumreclr(i,j,k) = sumreclr(i,j,k) +   &
                              reflectanceclr(i,j,k)*wtfac_p(i,j)
                      else if (daylight(i,j) ) then
                        sumtrclr(i,j,k) = sumtrclr(i,j,k) +   &
                              transmittance(i,j,k)*wtfac_p(i,j)
                        sumtr_dir_clr(i,j,k) =    &
                           sumtr_dir_clr(i,j,k) + tr_dir(i,j,k)*&
                                                          wtfac_p(i,j)
                        sumreclr(i,j,k) = sumreclr(i,j,k) +   &
                                      reflectance(i,j,k)*wtfac_p(i,j)
                      endif
                    end do
                  end do
                end do
              endif
            end do    ! end of gaussian loop
          end do  ! end of frequency points in the band loop
 
!----------------------------------------------------------------------
!    normalize the solar flux in the band to the appropriate value for  
!    the given total solar insolation.                                 
!---------------------------------------------------------------------
          solarflux_p(:,:) = fracday(:,:)*solflxband(nband)*ssolar/solflxtotal_local
 
          if (nband == visible_band_indx) then
            Sw_output%bdy_flx(:,:,1) =   &
                Sw_output%bdy_flx(:,:,1) + sumre(:,:,1)*   &
                                                      solarflux_p(:,:)
            Sw_output%bdy_flx(:,:,3) =    &
                Sw_output%bdy_flx(:,:,3) + sumtr(:,:,KERAD+1)*&
                                              solarflux_p(:,:) -  &
                                              sumre(:,:,KERAD+1)* &
                                              solarflux_p(:,:) 
          endif
          if (nband == onepsix_band_indx) then
             Sw_output%bdy_flx(:,:,2) =     &
                 Sw_output%bdy_flx(:,:,2) + sumre(:,:,1)*  &
                                                      solarflux_p(:,:)
             Sw_output%bdy_flx(:,:,4) =    &
                 Sw_output%bdy_flx(:,:,4) + sumtr(:,:,KERAD+1)*&
                                               solarflux_p(:,:) - &
                                               sumre(:,:,KERAD+1)*&
                                               solarflux_p(:,:)
          endif
          
!-------------------------------------------------------------------
!    calculate the band fluxes and heating rates.                       
!--------------------------------------------------------------------
          if (do_esfsw_band_diagnostics) then
            do k = KSRAD,KERAD+1
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  dfswband(i,j,k) = sumtr(i,j,k)* solarflux_p(i,j) 
                  ufswband(i,j,k) = sumre(i,j,k)* solarflux_p(i,j)
                end do
              end do
            end do
          endif
 
!----------------------------------------------------------------------
!    sum the band fluxes and heating rates to calculate the total       
!    spectral values.                                                  
!------------------------------------------------------------------
          do k = KSRAD,KERAD+1
            do j=JSRAD,JERAD
              do i=ISRAD,IERAD
                if (daylight(i,j) ) then
                  Sw_output%dfsw (i,j,k) =   &
                     Sw_output%dfsw(i,j,k) + sumtr(i,j,k)*&
                                                     solarflux_p(i,j)
                  Sw_output%ufsw (i,j,k) =   &
                     Sw_output%ufsw(i,j,k) + sumre(i,j,k)*  &
                                                     solarflux_p(i,j)
                  fswband(i,j,k) = ((sumre(i,j,k)*  &
                                        solarflux_p(i,j)) - &
                                       (sumtr(i,j,k)*    &
                                                    solarflux_p(i,j)))
                  Sw_output%fsw(i,j,k) = Sw_output%fsw(i,j,k) +&
                                                    fswband(i,j,k)
                endif
              end do
            end do
          end do
 
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              if (daylight(i,j) ) then
                Sw_output%dfsw_dir_sfc(i,j) =   &
                          Sw_output%dfsw_dir_sfc(i,j) +   &
                            sumtr_dir(i,j,KERAD+1)*solarflux_p(i,j)
                Sw_output%ufsw_dir_sfc(i,j) =   &
                          Sw_output%ufsw_dir_sfc(i,j) +   &
                            sumtr_dir_up(i,j)*solarflux_p(i,j)
                Sw_output%ufsw_dif_sfc(i,j) =   &
                           Sw_output%ufsw_dif_sfc(i,j) +   &
                               sumre(i,j,KERAD+1)*solarflux_p(i,j)
              endif
            end do
          end do

          if (nband > NIRBANDS) then
            do j=JSRAD,JERAD
              do i=ISRAD,IERAD
                if (daylight(i,j) ) then
                  Sw_output%dfsw_vis_sfc(i,j) =   &
                        Sw_output%dfsw_vis_sfc(i,j) +   &
                              sumtr(i,j,KERAD+1)*solarflux_p(i,j)
                  Sw_output%ufsw_vis_sfc(i,j) =   &
                         Sw_output%ufsw_vis_sfc(i,j) +   &
                               sumre(i,j,KERAD+1)*solarflux_p(i,j)
                  Sw_output%dfsw_vis_sfc_dir(i,j) =   &
                          Sw_output%dfsw_vis_sfc_dir(i,j) +   &
                            sumtr_dir(i,j,KERAD+1)*solarflux_p(i,j)
                  Sw_output%ufsw_vis_sfc_dir(i,j) =   &
                          Sw_output%ufsw_vis_sfc_dir(i,j) +   &
                            sumtr_dir_up(i,j)*solarflux_p(i,j)
                  Sw_output%ufsw_vis_sfc_dif(i,j) =   &
                           Sw_output%ufsw_vis_sfc_dif(i,j) +   &
                               sumre(i,j,KERAD+1)*solarflux_p(i,j)
                endif
              end do
            end do
          endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
          do k = KSRAD,KERAD
            do j=JSRAD,JERAD
              do i=ISRAD,IERAD
                if (daylight(i,j) ) then
                  hswband(i,j,k) = (fswband(i,j,k+1) -    &
                                    fswband(i,j,k) )*gocpdp(i,j,k)
                  Sw_output%hsw(i,j,k) =    &
                        Sw_output%hsw(i,j,k) + hswband(i,j,k)
                endif
              end do
            end do
          end do

!----------------------------------------------------------------------
!    calculate the band fluxes and heating rates.                       
!---------------------------------------------------------------------
          if (nprofile == 1) then  ! clr sky need be done only for 
                                     ! first cloud profile
            if (do_totcld_forcing) then
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  if (daylight(i,j) ) then
                    Sw_output%dfsw_dir_sfc_clr(i,j) =   &
                       Sw_output%dfsw_dir_sfc_clr(i,j) +   &
                        sumtr_dir_clr(i,j,KERAD+1)*solarflux_p(i,j)
                  endif
                end do
              end do
              if (nband > NIRBANDS) then
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    if (daylight(i,j) ) then
                      Sw_output%dfsw_vis_sfc_clr(i,j) =   &
                         Sw_output%dfsw_vis_sfc_clr(i,j) +   &
                            sumtrclr(i,j,KERAD+1)*solarflux_p(i,j)
                    endif
                  end do
                end do
              endif
              if (nband == visible_band_indx) then
                Sw_output%bdy_flx_clr(:,:,1) =      &
                         sumreclr(:,:,1)*solarflux_p(:,:)
                Sw_output%bdy_flx_clr(:,:,3) =    &
                        sumtrclr(:,:,KERAD+1)*solarflux_p(:,:) - &
                        sumreclr(:,:,KERAD+1)*solarflux_p(:,:) 
              endif
              if (nband == onepsix_band_indx) then
                Sw_output%bdy_flx_clr(:,:,2) =    &
                        sumreclr(:,:,1)*solarflux_p(:,:)
                Sw_output%bdy_flx_clr(:,:,4) =    &
                       sumtrclr(:,:,KERAD+1)*solarflux_p(:,:)  -  &
                       sumreclr(:,:,KERAD+1)*solarflux_p(:,:) 
              endif
          
              if (do_esfsw_band_diagnostics) then
                do k = KSRAD,KERAD+1
                  do j=JSRAD,JERAD
                    do i=ISRAD,IERAD
                      dfswbandclr(i,j,k) =     &
                                 sumtrclr(i,j,k)*solarflux_p(i,j)
                      ufswbandclr(i,j,k) =    &
                                 sumreclr(i,j,k)*solarflux_p(i,j)
                    end do
                  end do
                end do
              endif

!----------------------------------------------------------------------c
!    sum the band fluxes and heating rates to calculate the total     
!    spectral values.                                                 
!----------------------------------------------------------------------c
              do k = KSRAD,KERAD+1
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    Sw_output%dfswcf(i,j,k) =    &
                            Sw_output%dfswcf(i,j,k) +  &
                                  sumtrclr(i,j,k)*solarflux_p(i,j)
                    Sw_output%ufswcf(i,j,k) =      &
                            Sw_output%ufswcf(i,j,k) +  &
                                  sumreclr(i,j,k)*solarflux_p(i,j)
                    fswbandclr(i,j,k) =    &
                          (sumreclr(i,j,k)*solarflux_p(i,j)) - &
                          (sumtrclr(i,j,k)*solarflux_p(i,j))
                    Sw_output%fswcf(i,j,k) =    &
                                  Sw_output%fswcf(i,j,k) +    &
                                                fswbandclr(i,j,k)
                  end do
                end do
              end do

!----------------------------------------------------------------------c
!    sum the band fluxes and heating rates to calculate the total    
!    spectral values.                                               
!----------------------------------------------------------------------c
              do k = KSRAD,KERAD
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                     hswbandclr(i,j,k) =    &
                              (fswbandclr(i,j,k+1) -      &
                                  fswbandclr(i,j,k))*gocpdp(i,j,k)
                    Sw_output%hswcf(i,j,k) =   &
                                 Sw_output%hswcf(i,j,k) +    &
                                                 hswbandclr(i,j,k)
                  end do
                end do
              end do
            endif
          endif ! (nprofile == 1)
        end do   ! (end of band loop)
      end do   ! (end of nprofile loop)

!----------------------------------------------------------------------
!    if the ica calculation was being done, the fluxes and heating rates
!    which have been summed over nprofiles cloud profiles must be 
!    averaged.
!------------------------------------------------------------------
      if (do_ica_calcs) then
        do j=JSRAD,JERAD
          do i=ISRAD,IERAD
            Sw_output%dfsw_dir_sfc (i,j) =   &
                      Sw_output%dfsw_dir_sfc(i,j)*profiles_inverse
            Sw_output%ufsw_dif_sfc (i,j) =   &
                      Sw_output%ufsw_dif_sfc(i,j)*profiles_inverse
            Sw_output%dfsw_vis_sfc (i,j) =   &
                      Sw_output%dfsw_vis_sfc(i,j)*profiles_inverse
            Sw_output%ufsw_vis_sfc (i,j) =   &
                      Sw_output%ufsw_vis_sfc(i,j)*profiles_inverse
            Sw_output%dfsw_vis_sfc_dir (i,j) =   &
                      Sw_output%dfsw_vis_sfc_dir(i,j)*profiles_inverse
            Sw_output%ufsw_vis_sfc_dif (i,j) =   &
                      Sw_output%ufsw_vis_sfc_dif(i,j)*profiles_inverse
            Sw_output%bdy_flx(i,j,:) =  &
                      Sw_output%bdy_flx(i,j,:)*profiles_inverse
          end do
        end do
        do k = KSRAD,KERAD+1
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              Sw_output%dfsw (i,j,k) = Sw_output%dfsw(i,j,k)*  &
                                       profiles_inverse
              Sw_output%ufsw (i,j,k) = Sw_output%ufsw(i,j,k)*  &
                                       profiles_inverse
              Sw_output%fsw(i,j,k) = Sw_output%fsw(i,j,k)*  &
                                     profiles_inverse
            end do
          end do
        end do
        do k = KSRAD,KERAD
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              Sw_output%hsw(i,j,k) = Sw_output%hsw(i,j,k)*  &
                                       profiles_inverse
            end do
          end do
        end do
      endif

      do j=JSRAD,JERAD
        do i=ISRAD,IERAD
          if (daylight(i,j) ) then
            Sw_output%dfsw_dif_sfc(i,j ) =   &
                              Sw_output%dfsw(i,j,KERAD+1 ) -   &
                                      Sw_output%dfsw_dir_sfc(i,j)
          endif
        end do
      end do

      if (do_totcld_forcing) then
        do j=JSRAD,JERAD
          do i=ISRAD,IERAD
            if (daylight(i,j) ) then
              Sw_output%dfsw_dif_sfc_clr(i,j) =   &
                                Sw_output%dfswcf(i,j,KERAD+1) -   &
                                Sw_output%dfsw_dir_sfc_clr(i,j)
            endif
          end do
        end do
      endif

      do j=JSRAD,JERAD
        do i=ISRAD,IERAD
          if (daylight(i,j) ) then
            Sw_output%dfsw_vis_sfc_dif(i,j) =   &
                                Sw_output%dfsw_vis_sfc(i,j) -   &
                                Sw_output%dfsw_vis_sfc_dir(i,j) 
          endif
        end do
      end do

!--------------------------------------------------------------------
!    convert sw fluxes to cgs and then back to  mks units.
!---------------------------------------------------------------------
      Sw_output%fsw(:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%fsw(:,:,:))
      Sw_output%dfsw(:,:,:) =    &
                            1.0E-03*(1.0E+03*Sw_output%dfsw(:,:,:))
      Sw_output%ufsw(:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%ufsw(:,:,:))
      if (do_totcld_forcing) then
        Sw_output%fswcf(:,:,:) =   &
                            1.0E-03*(1.0E+03*Sw_output%fswcf(:,:,:))
        Sw_output%dfswcf(:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%dfswcf(:,:,:))
        Sw_output%ufswcf(:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%ufswcf(:,:,:))
      endif

!---------------------------------------------------------------------


end subroutine swresf



!####################################################################

subroutine esfsw_driver_end

!---------------------------------------------------------------------
!    esfsw_driver_end is the destructor for esfsw_driver_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('esfsw_driver_mod',   &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    close out the modules that this module initialized.
!--------------------------------------------------------------------
      call esfsw_utilities_end
      call esfsw_parameters_end

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

end subroutine esfsw_driver_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#################################################################
! <SUBROUTINE NAME="compute_gas_props">
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
!   call compute_gas_props (press, pflux, rh2o, deltaz,
!                           qo3, rrvco2, rrvch4, rrvn2o, cosz,
!                           daylight, gasopdep)
!  </TEMPLATE>
!  <IN NAME="press" TYPE="real">
!    layer pressures (mks units) (does not include sfc)
!  </IN>
!  <IN NAME="pflux" TYPE="real">
!    interface pressures (mks units)
!  </IN>
!  <IN NAME="rh2o" TYPE="real">
!    h2o mixing ratio
!  </IN>
!  <IN NAME="deltaz" TYPE="real">
!    layer thickness in meters
!  </IN>
!  <IN NAME="cosz" TYPE="real">
!    cosine of the zenith angle
!  </IN>
! </SUBROUTINE>

   subroutine compute_gas_props (press, pflux, temp, rh2o, deltaz, &
                                 qo3, rrvco2, rrvch4, rrvn2o, &
                                 cosz, daylight, gasopdep)

!----------------------------------------------------------------------
!    comput uses the delta-eddington technique in conjunction with a    
!    multiple-band parameterization for h2o+co2+o2+o3 absorption to   
!    derive solar fluxes and heating rates.                             
!    notes: drops are assumed if temp>273.15K, ice crystals otherwise.
!-------------------------------------------------------------------

real, dimension(:,:,:),        intent(in)    :: press, pflux, temp, rh2o, deltaz, qo3
real, dimension(:,:),          intent(in)    :: cosz
real,                          intent(in)    :: rrvco2, rrvch4, rrvn2o
logical, dimension(:,:),       intent(in)    :: daylight
real, dimension(:,:,:,:,:),    intent(out)   :: gasopdep              


!-------------------------------------------------------------------
!  intent(in) variables:
!
!      press          layer pressures (mks units) (does not include sfc)
!      pflux          interface pressures (mks units)
!      rh2o           h2o mixing ratio
!      deltaz         layer thickness in meters
!      rrvco2         carbon dioxide
!      rrvch4         methane
!      rrvn2o         nitrous oxide
!      cosz           cosine of the zenith angle
!                                                                 
!   intent(inout) variables:
!
!      gasopdep
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables:
 

      real, dimension (size(press,3))  :: &
                   efftauo2,   efftauco2,   efftauch4,   efftaun2o, &
                   wh2ostr,    wo3,         wo2,         quenchfac, &
                   opdep,      delpdig,     deltap,      tco2,    &
                   tch4,       tn2o,        to2,         wh2o,    &
                   wh2octms,   wh2octmf     

           
      real, dimension (size(pflux,3))  :: &
            alphaco2,        alphaco2str,    alphao2,          &
            alphao2str,      alphach4,       alphach4str,      &
            alphan2o,        alphan2ostr,    scale,            &
            scalestr,        totco2,         totco2str,        &
            toto2,           toto2str,       totch4,           &
            totch4str,       totn2o,         totn2ostr,        &
            pflux_mks,       z

      real, dimension (size(press,3)) :: sumwco2, sumwco2str
      real, dimension (size(pflux,3)) :: dnfitco2, sumdnfitco2
      real, dimension (size(press,3)) :: tauco2

      real :: cosangsolar
      real :: denom
      real :: wtquench
      real :: h2o_conv

      integer  :: j, i, k, ng, nband, kq
      integer  :: np, nf
      integer  :: israd, jsrad, ierad, jerad, ksrad, kerad


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
      ierad = size(press,1)
      jerad = size(press,2)
      kerad = size(press,3)

!---------------------------------------------------------------------
!    initialize local variables.                                        
!------------------------------------------------------------------
      alphaco2   (1) = 0.0
      alphaco2str(1) = 0.0
      alphao2    (1) = 0.0
      alphao2str (1) = 0.0
      alphach4   (1) = 0.0
      alphach4str(1) = 0.0
      alphan2o   (1) = 0.0
      alphan2ostr(1) = 0.0
      efftauo2(:) = 0.0
      efftauCo2(:) = 0.0
      efftaun2o(:) = 0.0
      efftauch4(:) = 0.0
      h2o_conv = 3.30E16
      wh2octms(:) = 0 
      wh2octmf(:) = 0

      do j = JSRAD,JERAD
        do i = ISRAD,IERAD
          if ( daylight(i,j) ) then

!----------------------------------------------------------------------c
!    define pressure related quantities, pressure is in mks units. 
!----------------------------------------------------------------------c

            do k = KSRAD+1,KERAD+1
              deltap(k-1) = pflux(i,j,k) - pflux(i,j,k-1)
              delpdig(k-1) = deltap(k-1)/ GRAV
              scalestr(k) = pflux(i,j,k) 
              scale(k) = scalestr(k)*pflux(i,j,k)/pstd_mks
            end do
 

            do k = KSRAD,KERAD
              wh2ostr(k) = rh2o(i,j,k)*delpdig(k)
              wo3(k)     = qo3(i,j,k)*delpdig(k)
              wo2(k) = o2mixrat*(WTMO2/WTMAIR)*delpdig(k)
            end do
  
!---------------------------------------------------------------------
!    if quenching factor effects are desired, calculate the height above
!    the surface of the model flux levels.
!---------------------------------------------------------------------
            if (do_quench) then
              z(KERAD+1) = 0.0
              do k = KERAD,KSRAD,-1
                z(k) = z(k+1) + deltaz(i,j,k)
              end do
          
!---------------------------------------------------------------------
!    define the quenching factor for each grid point.
!---------------------------------------------------------------------
              do k = KSRAD,KERAD
                if (z(k) < co2_quenchfac_height(1) ) then
                  quenchfac(k) = 1.0
                else if (z(k) > co2_quenchfac_height(30) ) then 
                  quenchfac(k) = 0.037
                else
                  do kq = 1,29
                    if (z(k) > co2_quenchfac_height(kq) .and. &
                        z(k) <= co2_quenchfac_height(kq+1)) then
                      wtquench = (z(k) - co2_quenchfac_height(kq))/ &
                                 (co2_quenchfac_height(kq+1) - &
                                  co2_quenchfac_height(kq))
                      quenchfac(k) = (1. - wtquench)*   &
                                           co2_quenchfac(kq) +   &
                                      wtquench*co2_quenchfac(kq+1)
                      exit
                    endif
                  end do
                endif
              end do
            else
              quenchfac(:) = 1.0
            endif !(do_quench)


            do ng = 1,NSOLWG
              cosangsolar = cosz(i,j)
              if (cosangsolar == 0.0) cosangsolar = 1.0

!----------------------------------------------------------------------c
!    define the scaled and unscaled co2 and o2 pathlengths in 
!    centimeter-atm, and the unscaled h2o and o3 amounts in   
!    kgrams/meter**2. 
!    cm-atm needed as units because of c2co2 having those units.
!----------------------------------------------------------------------c
              denom = 1.0/(GRAV*rhoair*cosangsolar*2.0)
              do k = KSRAD+1,KERAD+1
                totco2(k) = 1.0E+02*rrvco2*scale(k)*denom       
                totco2str(k) = 2.0E+02*rrvco2*scalestr(k)*denom     
                toto2(k) = 1.0E+02*o2mixrat*scale(k)*denom      
                toto2str(k) = 2.0E+02*o2mixrat*scalestr(k)*denom      
                if (do_ch4_sw_effects) then
                  totch4(k) = 1.0E+02*rrvch4*scale(k)*denom     
                  totch4str(k) = 2.0E+02*rrvch4*scalestr(k)*denom     
                endif
                if (do_n2o_sw_effects) then
                  totn2o(k) = 1.0E+02*rrvn2o*scale(k)*denom     
                  totn2ostr(k) = 2.0E+02*rrvn2o*scalestr(k)*denom      
                endif
              end do

              np = 0
              do nband = 1, NBANDS

!-------------------------------------------------------------------
!    define the h2o scaled gas amounts in kgrams/meter**2
!    
!    Define the h2o self and foriegn continuum scaled gas amounts
!    in molecules/(cm**2*atm) 
!    See equation 5b of Paynter et al. 2009 
!    Optical depth self continuum = Number of moleclues in path*
!    (vapor pressure/101325)*Continuum Coefficent*
!    temperature dependence of continuum cofficent
!    Optical depth foriegn continuum = Number of moleclues in path*
!    (dry atm pressure/101325)*Continuum Coefficent*
!    temperature dependence of contiuum cofficent 
!
!    We firstly have to convert from mass density
!    (i.e. rh2o*delpdig) to number density by mutiplying by 
!    (R_water/K_b) and then mulitply this value by vapor 
!     pressure of H2O/1atm for self continuum and pressure of 
!     AIR/1atm for foreign continuum.
!     We will assume presssure of H2O = sphum*Pressure(atm)/0.622
!     and use 1E-4 converts from m2 to cm2  
!    
!    For self continuum: 
!    Number of moleclues in path*(vapor pressure/101325) = 
!    press*sphum*(deltaP/g)*(R_water/K_b)*(1E-4*(sphum/0.622)/101325) 
!    Define h2o_conv as  (R_water/K_b)*1E-4/101325 = 3.30E16 
!    sphum**2*press/0.622*(deltaP/g)*h2o_conv  
!    
!    For foreign continuum:  
!    sphum*press*(deltaP/g)*h2o_conv  
!                            
!    The temperature dependence of the self continuum has the form:
!    exp(TD*(1/T-1/296)) 
!    We have assumed a value of TD of 1500, based up Paynter 
!    and ramaswamy 2014 JGR. 
!    The temperature dependence of the foriegn continuum has the form
!    296/T). 
!     
!---------------------------------------------------------------------
                if (nband <= nh2obands) then
     
                  do k = KSRAD,KERAD
                    wh2o(k) = rh2o(i,j,k)*delpdig(k)*   &
                        exp(powph2o(nband)*alog(press(i,j,k)*p0h2o(nband)))
                    if (do_sw_continuum) then 
                       wh2octms(k) = (h2o_conv/0.622)*press(i,j,k)*rh2o(i,j,k)**2 &
                            *delpdig(k)*exp(1500.*((1./temp(i,j,k)) -3.38E-3))&
                            *(296./temp(i,j,k))
                       wh2octmf(k) = h2o_conv*press(i,j,k)*rh2o(i,j,k)&
                            *delpdig(k)*(296./temp(i,j,k)) - wh2octms(k)     
                    end if 
                  end do

!---------------------------------------------------------------------
!    calculate the "effective" co2, o2, ch4 and n2o gas optical depths 
!    for the appropriate absorbing bands.                               
!    note: for large zenith angles, alpha can exceed 1. In this case,a
!    the optical depths are set to the previous layer values.          
!-------------------------------------------------------------------
                
                  if ((nbands == 18 .and. nfrqpts == 38) .or. reproduce_ulm) then
                    if ( c1co2(nband).ne.1.0E-99 ) then
                      do k = KSRAD+1,KERAD+1
                        if (totco2(k) < totco2max(nband) .and.  &
                            totco2str(k) < totco2strmax(nband))  then
                          alphaco2(k) =     &
                               c1co2(nband)*exp(c3co2(nband)* &
                                   alog((totco2(k) + c2co2(nband))))  -  &
                                                            c4co2(nband)
                          alphaco2str(k) = &
                            c1co2str(nband)*exp(c3co2str(nband)*  &
                              alog((totco2str(k) + c2co2str(nband)))) - &
                                                          c4co2str(nband)
                          tco2(k-1) =      &
                               (1.0 - alphaco2(k))*   &
                                              (1.0 - alphaco2str(k))/ &
                               ((1.0 - alphaco2(k-1))*    &
                                              (1.0 - alphaco2str(k-1)))
                          efftauco2(k-1) = -cosangsolar*alog( tco2(k-1))
                        else if (k > KSRAD+1) then
                          efftauco2(k-1) = efftauco2(k-2)
                        else
                          efftauco2(k-1) = 0.0
                          endif
                      end do
                    else    !( c1co2(nband).ne.1.0E-99 ) 
                      efftauco2(:) = 0.0
                    end if  !( c1co2(nband).ne.1.0E-99 ) 

                  else  ! nbands

                    if (do_co2_sw_effects) then
                      if (nband.le.NCO2BANDS) then
                        do k = KSRAD,KERAD
                          if (k.eq.1) then
                            sumwco2(k) = 0.1 * rrvco2 * (WTMCO2/WTMAIR) * delpdig(k) *   &
                                         exp( powpco2(nband) * alog( press(i,j,k)*p0co2(nband) ) )
                            sumwco2str(k) = 0.1 * rrvco2*(WTMCO2/WTMAIR)*delpdig(k)
                          else
                            sumwco2(k) = sumwco2(k-1) + &
                                         0.1 * rrvco2 * (WTMCO2/WTMAIR) * delpdig(k) * &
                                         exp( powpco2(nband) * alog( press(i,j,k)*p0co2(nband) ) )
                            sumwco2str(k) = sumwco2str(k-1) + &
                                            0.1 * rrvco2 * (WTMCO2/WTMAIR) * delpdig(k)
                          endif
                        end do
                        sumdnfitco2 = 0.0
                        do nf = 1,nfreqptsco2(nband)
                          dnfitco2(1) = 1.0
                          if (.not.strtermco2(nf,nband)) then
                            do k = KSRAD,KERAD
                              tauco2(k) = kco2(nf,nband) * sumwco2(k)
                            end do
                          endif
                          if (strtermco2(nf,nband)) then
                            do k = KSRAD,KERAD
                              tauco2(k) = kco2(nf,nband) * sumwco2str(k)
                            end do
                          endif
                          do k = KSRAD,KERAD
                            dnfitco2(k+1) = exp(-tauco2(k)/cosangsolar)
                          end do

                          sumdnfitco2(1) = 1.0
                          do k = KSRAD+1,KERAD+1
                            sumdnfitco2(k) = sumdnfitco2(k) + wtfreqco2(nf,nband) * dnfitco2(k)
                          end do

                        end do ! nf

                        do k = KSRAD,KERAD
                          efftauco2(k) = -cosangsolar * alog( sumdnfitco2(k+1)/sumdnfitco2(k) )
                        end do
 
                      else  ! nband.le.NCO2BANDS
                        efftauco2(:) = 0.0
                      endif

                    else ! do_co2_sw_effects
                      efftauco2(:) = 0.0
                    endif
                  endif ! nband

                  if (do_ch4_sw_effects) then
                    if (c1ch4(nband).ne.1.0E-99 ) then
                      do k = KSRAD+1,KERAD+1
                        if (totch4(k) < totch4max(nband) .and.  &
                              totch4str(k) < totch4strmax(nband))  then
                           alphach4(k) =    &
                              c1ch4(nband)*exp(c3ch4(nband)*&
                              alog((totch4(k) + c2ch4(nband))))  -   &
                                                           c4ch4(nband)
                           alphach4str(k) = &
                            c1ch4str(nband)*exp(c3ch4str(nband)*  &
                            alog((totch4str(k) + c2ch4str(nband)))) - &
                                                        c4ch4str(nband)
                           tch4(k-1) = &
                                  (1.0 - alphach4(k))*    &
                                           (1.0 - alphach4str(k))/ &
                                   ((1.0 - alphach4(k-1))*   &
                                           (1.0 - alphach4str(k-1)))
                           efftauch4(k-1) = -cosangsolar*alog(tch4(k-1))
                        else if (k > KSRAD+1) then
                           efftauch4(k-1) = efftauch4(k-2)
                        else
                           efftauch4(k-1) = 0.0
                        end if
                      end do
                    else    !( c1ch4(nband).ne.1.0E-99 )
                      efftauch4(:) = 0.0
                    end if  !( c1ch4(nband).ne.1.0E-99 )
                  else    !do_ch4 = .false.
                    efftauch4(:) = 0.0
                  end if

                  if (do_n2o_sw_effects) then
                    if ( c1n2o(nband).ne.1.0E-99 ) then
                      do k = KSRAD+1,KERAD+1
                        if (totn2o(k) < totn2omax(nband) .and.  &
                            totn2ostr(k) < totn2ostrmax(nband)) then
                          alphan2o(k) = &
                               c1n2o(nband)*exp(c3n2o(nband)* &
                                  alog((totn2o(k) +c2n2o(nband)))) -  &
                                                         c4n2o(nband)
                          alphan2ostr(k) = &
                            c1n2ostr(nband)*exp(c3n2ostr(nband)*  &
                            alog((totn2ostr(k) + c2n2ostr(nband)))) -  &
                                                 c4n2ostr(nband)
                          tn2o(k-1) = &
                                    (1.0 - alphan2o(k)) *  &
                                    (1.0 - alphan2ostr(k))/ &
                                  (( 1.0 - alphan2o(k-1)) *  &
                                    (1.0 - alphan2ostr(k-1)))
                          efftaun2o(k-1) = -cosangsolar*alog(tn2o(k-1))
                        else if (k > KSRAD+1) then
                          efftaun2o(k-1) = efftaun2o(k-2)
                        else
                          efftaun2o(k-1) = 0.0
                        end if
                      end do
                    else    !( c1n2o(nband).ne.1.0E-99 )
                      efftaun2o(:) = 0.0
                    end if  !( c1n2o(nband).ne.1.0E-99 )
                  else  !do_n2o = .false.
                    efftaun2o(:) = 0.0
                  end if



                  if (reproduce_ulm) then
                    if (do_o2_sw_effects .and. c1o2(nband).ne.1.0E-99 ) then
                      do k = KSRAD+1,KERAD+1                     
                        if (toto2(k) .lt. toto2max(nband) .and.   &
                            toto2str(k) .lt. toto2strmax(nband)) then
                          alphao2(k) = c1o2(nband)*exp( c3o2(nband)* &
                                       alog((toto2(k) + c2o2(nband)))) - &
                                                        c4o2(nband)
                          alphao2str(k) = &
                                  c1o2str(nband)*exp(c3o2str(nband)*  &
                                      alog((toto2str(k) + c2o2str(nband)))) &
                                                              - c4o2str(nband)
                          to2(k-1) = &
                                     (1.0 - alphao2(k))*  &
                                     (1.0 - alphao2str(k) )/ &
                                    ((1.0 - alphao2(k-1)) *  &
                                     (1.0 - alphao2str(k-1)))

                          efftauo2(k-1) = -cosangsolar*alog(to2(k-1))

                        else if (k.gt.KSRAD+1) then
                          !write(3322,*) '1S', nband, toto2str(k), alphao2str(k)
                          !write(3322,*) '1', nband, toto2(k), alphao2(k),efftauo2(k-2)
                          efftauo2(k-1) = efftauo2(k-2)
                        else
                          !write(3322,*) '2S', toto2str(k), alphao2str(k)
                          !write(3322,*) '2', toto2(k), alphao2(k),efftauo2(k-2)
                          efftauo2(k-1) = 0.0 
                        end if 
                      end do

                    else   !  ( c1o2(nband).ne.1.0E-99 ) 
                      efftauo2(:) = 0.0
                    end if  !  ( c1o2(nband).ne.1.0E-99 ) 
                  

                  else ! reproduce_ulm
                  

                    if (do_o2_sw_effects .and. c1o2(nband).ne.1.0E-99  .and. c1o2str(nband).ne.1.0E-99) then
   
                      do k = KSRAD+1,KERAD+1
                     
                        if (toto2(k) .lt. toto2max(nband) .and.   &
                            toto2str(k) .lt. toto2strmax(nband)) then
                          alphao2(k) = c1o2(nband)*exp( c3o2(nband)* &
                                       alog((toto2(k) + c2o2(nband)))) - &
                                                        c4o2(nband)
                          alphao2str(k) = &
                                  c1o2str(nband)*exp(c3o2str(nband)*  &
                                      alog((toto2str(k) + c2o2str(nband)))) &
                                                              - c4o2str(nband)
                          to2(k-1) = &
                                     (1.0 - alphao2(k))*  &
                                     (1.0 - alphao2str(k) )/ &
                                    ((1.0 - alphao2(k-1)) *  &
                                     (1.0 - alphao2str(k-1)))

                          efftauo2(k-1) = -cosangsolar*alog(to2(k-1))

                        else if (k.gt.KSRAD+1) then
                          efftauo2(k-1) = efftauo2(k-2)
                        else
                          efftauo2(k-1) = 0.0 
                        end if 

                      end do

                    else  if (do_o2_sw_effects .and. c1o2(nband).ne.1.0E-99 .and. c1o2str(nband).eq.1.0E-99) then
 
                      do k = KSRAD+1,KERAD+1
                     
                        if (toto2(k) .lt. toto2max(nband)) then
                          alphao2(k) = c1o2(nband)*exp( c3o2(nband)* &
                                       alog((toto2(k) + c2o2(nband)))) - &
                                                        c4o2(nband)
                
                          to2(k-1) = &
                                     (1.0 - alphao2(k))*  &
                                          (1.0)/ &
                                      ((1.0 - alphao2(k-1)) *  &
                                                  (1.0))

                          efftauo2(k-1) = -cosangsolar*alog(to2(k-1))

                        else if (k.gt.KSRAD+1) then
                          efftauo2(k-1) = efftauo2(k-2)
                        else
                          efftauo2(k-1) = 0.0 
                        end if 

                      end do
                                     
                    else   !  ( c1o2(nband).ne.1.0E-99 ) 
                      efftauo2(:) = 0.0
                    end if  !  ( c1o2(nband).ne.1.0E-99 ) 		   
   
                  end if ! reproduce_ulm


                else ! if (nband > nh2obands) then

                  if (reproduce_ulm) then
                  else
                    efftauo2(:) = 0.0
                    efftauCo2(:) = 0.0
                    efftaun2o(:) = 0.0
                    efftauch4(:) = 0.0
                  end if  
                end if  ! (nband <= nh2obands)
 
!---------------------------------------------------------------------
!    calculate the "effective" o2 gas optical depths for the Schuman- 
!    Runge band.                                                        
!-------------------------------------------------------------------
                if ( nband.EQ.NBANDS ) then
                  do k = KSRAD+1,KERAD+1
                    if ( toto2str(k).lt.toto2strmaxschrun) then
                      alphao2str(k) =  &
                           c1o2strschrun*exp( c3o2strschrun*&
                              alog((toto2str(k) + c2o2strschrun))) - &
                                                       c4o2strschrun
                      to2(k-1) = &
                          (1.0 - alphao2str(k))/(1.0 - alphao2str(k-1)) 
                      efftauo2(k-1) =  -cosangsolar*alog(to2(k-1) )
                      if (do_herzberg) then
                        efftauo2(k-1) = efftauo2(k-1) +     &
                                                 wo2(k-1)*herzberg_fac
                      endif
                    else if (k.gt.KSRAD+1) then
                      efftauo2(k-1) = efftauo2(k-2)
                    else
                      efftauo2(k-1) = 0.0
                    end if
                  end do
                end if

                do nf =1,nfreqpts(nband)
                  np = np + 1

!---------------------------------------------------------------------
!    define the h2o + o3 gas optical depths.                           
!--------------------------------------------------------------------
                  if (strterm(np)) then
                    opdep(:) = kh2o(np)*wh2ostr(:) + ko3(np)*wo3(:)
                  else
                    opdep(:) = kh2o(np)*wh2o(:) + ko3(np)*wo3(:)  
                  end if

                  if (do_sw_continuum) then 
                    gasopdep(i,j,:,np,ng) =    &
                          opdep(:) + quenchfac(:)*efftauco2(:) +   &
                          efftauo2(:) + efftauch4(:) + efftaun2o(:) &
                          +  wh2octms(:)*kctms(np) +wh2octmf(:)*kctmf(np)
                  else
                    gasopdep(i,j,:,np,ng) =    &
                          opdep(:) + quenchfac(:)*efftauco2(:) +   &
                          efftauo2(:) + efftauch4(:) + efftaun2o(:)
                  end if

                end do  ! (nf loop)
              end do   ! (nband loop)
            end do  ! (ng loop)
          endif  ! (daylight)
        end do ! (i loop)
      end do ! (j loop)

!---------------------------------------------------------------------


end subroutine compute_gas_props


!#####################################################################
!<SUBROUTINE NAME="adding">
! <OVERVIEW>
!  Subroutine that implements doubling and adding technique to combine
!  multiple atmospheric layers to calculate solar fluxes
! </OVERVIEW>
! <DESCRIPTION>
!  This subroutine implements the standard doubling and adding
!  technique to combine reflectance and transmittance of multiple 
!  atmospheric layers to compute solar flux and heating rate.
! </DESCRIPTION>
! <TEMPLATE>
!  call adding ( ix, jx, kx, &
!                rlayerdir, tlayerdir, rlayerdif, tlayerdif,  &
!                tlayerde, sfcalb, calc_flag, reflectance,   &
!                transmittance)
! </TEMPLATE>
! <IN NAME="ix" TYPE="integer">
!  ix is the current longitudinal index in the physics cell being
!  integrated.
! </IN>
! <IN NAME="jx" TYPE="integer">
!  jx is the current latitudinal index in the physics cell being
!  integrated.
! </IN>
! <IN NAME="kx" TYPE="integer">
!  ix is the current vertical index in the physics cell being
!  integrated.
! </IN>
! <IN NAME="rlayerdir" TYPE="real">
!  layer reflectivity to direct incident beam
! </IN>
! <IN NAME="tlayerdir" TYPE="real">
!  layer transmissivity to direct incident beam
! </IN>
! <IN NAME="rlayerdif" TYPE="real">
!  layer reflectivity to diffuse incident beam
! </IN>
! <IN NAME="tlayerdir" TYPE="real">
!  layer transmissivity to diffuse incident beam
! </IN>
! <IN NAME="tlayerde" TYPE="real">
!  layer diffuse transmissivity to direct incident beam
! </IN>
! <IN NAME="sfcalb" TYPE="real">
!  surface albedo
! </IN>
! <IN NAME="calcflag" TYPE="integer">
!  flag to indicate columns where adding is to be done
! </IN>
! <OUT NAME="reflectance" TYPE="real">
!  diffuse reflectance at a level
! </OUT>
! <OUT NAME="transmittance" TYPE="real">
!  diffuse transmittance at a level
! </OUT>
!</SUBROUTINE>
!

subroutine adding (ix, jx, kx, rlayerdir, tlayerdir, rlayerdif,   &
                   tlayerdif, tlayerde, sfcalb_dir, sfcalb_dif,  &
                   calc_flag, reflectance, transmittance, tr_dir)
 
!-------------------------------------------------------------------
!    adding calculates the reflection and transmission at flux levels 
!    from the direct and diffuse values of reflection and transmission
!    in the corresponding layers using the adding method.           
!    references:                                                        
!    bowen, m.m., and v. ramaswamy, effects of changes in radiatively
!        active species upon the lower stratospheric temperatures.,    
!        j. geophys. res., 18909-18921, 1994.                         
!--------------------------------------------------------------------

integer, intent(in)                    :: ix, jx, kx
real, dimension(:,:,:),   intent(in)   :: rlayerdir, rlayerdif, &
                                          tlayerdir, tlayerdif, & 
                                          tlayerde
real, dimension (:,:),    intent(in)   :: sfcalb_dir, sfcalb_dif
logical, dimension (:,:), intent(in)   :: calc_flag
real, dimension(:,:,:),   intent(out)  :: reflectance, transmittance, &
                                          tr_dir

!-------------------------------------------------------------------
!  intent(in) variables:
!
!    ix,jx,kx        dimensions of current physics window            
!    rlayerdir       layer reflectivity to a direct incident beam      
!    tlayerdir       layer transmissivity to a direct incident beam   
!    rlayerdif       layer reflectivity to a diffuse incident beam  
!    tlayerdif       layer transmissivity to a diffuse incident beam  
!    tlayerde        layer transmissivity (non-scattered) to the direct 
!                    incident beam                                 
!    sfcalb_dir      surface albedo, direct beam 
!    sfcalb_dif      surface albedo, diffuse beam
!    calc_flag       flag to indicate columns where adding is to be 
!                    done. calculations not done in "dark" columns and 
!                    on clr sky pass in columns without any clouds.
!
!  intent(out) variables:
!
!    reflectance     reflectance of the scattered radiation at a level 
!    transmittance   transmittance of the scattered radiation at a level
!    tr_dir
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables:
 
      real, dimension (lbound(rlayerdir,3):ubound(rlayerdir,3)+1) ::  &
                            raddupdif2, raddupdir2

      real, dimension (lbound(rlayerdir,3):ubound(rlayerdir,3)  ) ::  &
                                      radddowndif2,  tadddowndir2

      real :: dm1tl2, dm2tl2, rdm2tl2, dm32, dm3r2, dm3r1p2, alpp2, &
              raddupdif2p, raddupdir2p, tlevel2p, radddowndifm, &
              tadddowndirm
      integer     ::  k, j, i

!-------------------------------------------------------------------
!   local variables:
!
!      raddupdif2
!      raddupdir2
!      tlevel2
!      radddowndif2
!      tadddowndir2
!      tlayerdif2
!      tlayerdir2
!      rlayerdif2
!      rlayerdir2
!      tlayerde2
!      dm1tl2
!      dm2tl2
!      rdm2tl2
!      dm32
!      dm3r2
!      dm3r1p2
!      alpp2
!      raddupdif2
!      raddupdir2p
!      tlevel2p
!      radddowndifm
!      tadddowndirm
!      i,j,k
!
!--------------------------------------------------------------------

!----------------------------------------------------------------------c
!    initialization for the surface layer.                           
!----------------------------------------------------------------------c
      do j=1,jx        
        do i=1,ix
          if (calc_flag(i,j) ) then
 
!------------------------------------------------------------------ 
!    add the inhomogeneous layers upward from the surface to the top of
!    the atmosphere.                                                  
!    radiation incident from above for diffuse beam, reflection of  
!    direct beam and conversion to diffuse.                           
!--------------------------------------------------------------------
            raddupdif2p = sfcalb_dif(i,j)
            raddupdir2p = sfcalb_dir(i,j)
            do k = kx, 1,-1
              dm2tl2    = tlayerdif(i,j,k)/(1.0 - rlayerdif(i,j,k)*  &
                          raddupdif2p )
              rdm2tl2    = dm2tl2*raddupdif2p     
              raddupdif2(k) = rlayerdif(i,j,k) + tlayerdif(i,j,k)*   &
                              rdm2tl2    
              raddupdir2(k) = rlayerdir(i,j,k) + tlayerde(i,j,k)*   &
                              raddupdir2p* dm2tl2 +   &     
                              (tlayerdir(i,j,k) - tlayerde(i,j,k))*  &
                              rdm2tl2   
              raddupdir2p = raddupdir2(k)
              raddupdif2p = raddupdif2(k)
            end do
 
!---------------------------------------------------------------------
!    define the direct transmittance. add the inhomogeneous layers 
!    downward from the second layer to the surface. radiation incident
!    from below for diffuse beam, transmission of direct beam and 
!    conversion to diffuse.                             
!-------------------------------------------------------------------
 
!--------------------------------------------------------------------
!    initialization for the first scattering layer.                   
!-------------------------------------------------------------------
            tlevel2p         = tlayerde(i,j,1)
            radddowndifm    =  rlayerdif(i,j,1)
            tadddowndirm    =  tlayerdir(i,j,1)
            do k= 2,kx    
              dm1tl2 = tlayerdif(i,j,k)/(1.0 - rlayerdif(i,j,k)*  &
                       radddowndifm)
              radddowndif2(k) = rlayerdif(i,j,k) + radddowndifm* &
                                tlayerdif(i,j,k)*dm1tl2      
              tadddowndir2(k) = tlevel2p*(tlayerdir(i,j,k) + &
                                rlayerdir(i,j,k)*radddowndifm* &
                                dm1tl2) + (tadddowndirm -  &
                                tlevel2p)*dm1tl2           

!---------------------------------------------------------------------
!    add downward to calculate the resultant reflectances and           
!    transmittances at flux levels.                                    
!------------------------------------------------------------------
              dm32  = 1.0/(1.0 - raddupdif2(k)*radddowndifm)
              dm3r2 = dm32*radddowndifm      
              dm3r1p2 = 1.0 + raddupdif2(k)*dm3r2   
              alpp2 = (tadddowndirm - tlevel2p)*dm32   
              transmittance(i,j,k) = (tlevel2p*(1.0 + raddupdir2(k)* &
                                      dm3r2) + alpp2)
              tr_dir(i,j,k) = tlevel2p
              reflectance(i,j,k) = (tlevel2p*raddupdir2(k)*   &
                                    dm3r1p2 + alpp2*   &
                                    raddupdif2(k))
              tlevel2p = tlevel2p*tlayerde (i,j,k) 
              radddowndifm = radddowndif2(k)
              tadddowndirm = tadddowndir2(k)
            end do
!! CORRECT ???
!           dm32  = 1.0/(1.0 - sfcalb(i,j)*radddowndifm)
            dm32          = 1.0/(1.0 - sfcalb_dif(i,j)*   &
                               radddowndifm       )
            dm3r2 = dm32*radddowndifm       
!! CORRECT ???
!           dm3r1p2 = 1.0 + sfcalb(i,j)*dm3r2         
            dm3r1p2          = 1.0 + sfcalb_dif(i,j) * dm3r2
            alpp2 = (tadddowndirm - tlevel2p)*dm32          
            transmittance(i,j,kx+1) = (tlevel2p*(1.0 +   &
!! CORRECT ???
!                                      sfcalb(i,j)* &
!12-08-03:  CHANGE THIS TO _dir as per SMF  sfcalb_dif(i,j)* &
                                       sfcalb_dir(i,j)* &
                                       dm3r2) + alpp2)
            tr_dir(i,j,kx+1) = tlevel2p
            reflectance(i,j,kx+1) = (tlevel2p*  &
!! CORRECT ???
!                                   sfcalb(i,j)*   &
                                    sfcalb_dir(i,j)* &
                                     dm3r1p2 + alpp2* &
!! CORRECT ???
!                                sfcalb(i,j) )
                                    sfcalb_dif(i,j))  
            reflectance(i,j,1) = raddupdir2p         
            transmittance(i,j,1) = 1.0
            tr_dir(i,j,1) = 1.0
          endif
        end do
      end do

!------------------------------------------------------------------


end subroutine adding 


!####################################################################
! <SUBROUTINE NAME="deledd">
!  <OVERVIEW>
!   Subroutine that calculates reflectivity and transmissivity in a
!   scattering layer using delta-eddington method
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine takes layer optical depth, single scattering abledo,
!   and asymmetry parameter, using delta-eddington method, to calculate
!   direct/diffuse reflectivity/transmissivity to direct/diffuse incident
!   radiation. The approximation uses the strong forward scattering of
!   aerosol particles.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call deledd (ix, jx, kx,  &
!                taustr, omegastr, gstr, cosang, ng , daylight,  &
!                rlayerdir, tlayerdir, rlayerdif, tlayerdif,   &
!                tlayerde,  cloud)
!  </TEMPLATE>
!  <IN NAME="ix" TYPE="integer">
!  ix is the current longitudinal index in the physics cell being
!  integrated.
!  </IN>
!  <IN NAME="jx" TYPE="integer">
!   jx is the current latitudinal index in the physics cell being
!   integrated.
!  </IN>
!  <IN NAME="kx" TYPE="integer">
!   ix is the current vertical index in the physics cell being
!   integrated.
!  </IN>
!  <IN NAME="taustr" TYPE="real">
!   the scaled optical depth, true optical depth normalized using
!   delta-eddington approximation
!  </IN>
!  <IN NAME="omegastr" TYPE="real">
!   the scaled single-scattering albedo
!  </IN>
!  <IN NAME="gstr" TYPE="real">
!   the scaled asymmetry factor
!  </IN>
!  <IN NAME="cosang" TYPE="real">
!   cosine of the solar zenith angle
!  </IN>
!  <IN NAME="ng" TYPE="real">
!   the number of gaussian angles to compute the diurnally    
!   averaged solar radiation (=1 unless lswg = true)
!  </IN>
!  <IN NAME="cloud" TYPE="real">
!   flag for existence of a cloud (used only in 'ovc' mode)
!  </IN>
!  <OUT NAME="rlayerdir" TYPE="real">
!   layer reflectivity to direct incident beam
!  </OUT>
!  <OUT NAME="tlayerdir" TYPE="real">
!   layer transmissivity to direct incident beam
!  </OUT>
!  <OUT NAME="rlayerdif" TYPE="real">
!   layer reflectivity to diffuse incident beam
!  </OUT>
!  <OUT NAME="tlayerdir" TYPE="real">
!   layer transmissivity to diffuse incident beam
!  </OUT>
!  <OUT NAME="tlayerde" TYPE="real">
!   layer diffuse transmissivity to direct incident beam
!  </OUT>
! </SUBROUTINE>
!
subroutine deledd (ix, jx, kx, taustr, omegastr, gstr, cosang, ng, &
                   daylight, rlayerdir, tlayerdir, tlayerde,   &
                   rlayerdif, tlayerdif, cloud)
 
!---------------------------------------------------------------------- 
!    deledd calculates the reflection and transmission in the 
!    scattering layers using the delta-eddington method.         
!    references:                                                   
!      joseph, j.h., w. wiscombe, and j.a. weinman, the delta-eddington
!      approximation for radiative flux transfer.,j. atmos. sci.,33,  
!      2452-2459, 1976.                                              
!-------------------------------------------------------------------

integer,                   intent(in)              :: ix, jx, kx
real, dimension(:,:,:),    intent(inout)           :: taustr, omegastr
real, dimension(:,:,:),    intent(in)              :: gstr
real, dimension(:,:),    intent(in)                ::  cosang
integer,                   intent(in)              :: ng
logical, dimension(:,:),   intent(in)              :: daylight
real, dimension(:,:,:),    intent(out)             :: rlayerdir,   &
                                                      tlayerdir,   &
                                                      tlayerde
real, dimension(:,:,:),    intent(inout), optional :: rlayerdif,   &
                                                      tlayerdif
logical, dimension(:,:,:), intent(in), optional    :: cloud          

!----------------------------------------------------------------------
!  intent(in) variables:
!
!    ix,jx,kx
!    gstr        the scaled asymmetry factor                       
!    cosang      the cosine of the solar zenith angle    
!    ng          the number of gaussian angles to compute the diurnally 
!                averaged solar radiation (=1 unless lswg = true)       
!    daylight
!
!  intent(inout) variables:
!
!    taustr      the scaled extinction optical depth                    
!    omegastr    the scaled single-scattering albedo               
!
!  intent(out) variables:
!
!    rlayerdir   the layer reflectivity to a direct incident beam      
!    tlayerdir   the layer transmissivity to a direct incident beam   
!    rlayerdif   the layer reflectivity to a diffuse incident beam   
!    tlayerdif   the layer transmissivity to a diffuse incident beam
!    tlayerde    the layer transmissivity (non-scattered) to the direct 
!                incident beam                                       
!
! intent(in),optional:
!
!    cloud       flag for existence of a cloud (used only in 'ovc' mode)
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      real        :: qq(7), rr(5), ss(8), tt(8), ww(4)
      real        :: rsum, tsum
      real        :: onedi3 = 1.0/3.0           
      real        :: twodi3 = 2.0/3.0             
      integer     :: k, i, ns, j, nn, ntot

      integer, dimension(ix, jx, kx) :: cld_index

      real,    dimension(ix)                  ::   &
                                          gstr2, taustr2, omegastr2, &
                                           cosangzk2, rlayerdir2,    &
                                           tlayerde2, tlayerdir2, &
                                           sumr, sumt


!----------------------------------------------------------------------
!  local variables:
!
!      qq
!      rr
!      ss
!      tt
!      ww
!      rsum
!      tsum
!      alpha
!      onedi3
!      twodi3
!      i,j,k
!      ns
!      nn
!      ntot
!
!---------------------------------------------------------------------

      do k=1,kx         
        do j=1,jx         

!---------------------------------------------------------------------
!    overcast sky mode. note: in this mode, the delta-eddington method
!    is performed only for spatial points containing a cloud.   
!-------------------------------------------------------------------
          nn = 0
          if (present(cloud)) then
            do i=1,ix          
              if (cloud(i,j,k) ) then
                nn = nn + 1
                cld_index(i,j,k) = nn
                gstr2(nn) = gstr(i,j,k)
                taustr2(nn) = taustr(i,j,k)
                omegastr2(nn) = omegastr(i,j,k)
                cosangzk2(nn) = cosang(i,j)

!----------------------------------------------------------------------
!    note: the following are done to avoid the conservative scattering 
!    case, and to eliminate floating point errors in the exponential 
!    calculations, respectively.                      
!----------------------------------------------------------------------c
                if (omegastr2(nn) >= 1.0) omegastr2(nn) = 9.9999999E-01
                if (taustr2(nn) >= 1.0E+02) taustr2(nn) = 1.0E+02
              endif
            end do

!----------------------------------------------------------------------c
!    clear sky mode. note: in this mode, the delta-eddington method is 
!    performed for all spatial points.                 
!----------------------------------------------------------------------c
          else
            do i=1,ix         
              if (daylight(i,j) ) then
                nn = nn + 1
                cld_index(i,j,k) = nn
                gstr2(nn) = gstr(i,j,k)
                taustr2(nn) = taustr(i,j,k)
                omegastr2(nn) = omegastr(i,j,k)
                cosangzk2(nn) = cosang(i,j   )

!----------------------------------------------------------------------c
!    note: the following are done to avoid the conservative scattering  
!    case, and to eliminate floating point errors in the exponential 
!    calculations, respectively.                    
!----------------------------------------------------------------------c
                if (omegastr2(nn) >= 1.0) omegastr2(nn) = 9.9999999E-01
                if (taustr2(nn) >= 1.0E+02) taustr2(nn) = 1.0E+02
              endif
            end do
          endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          ntot = nn

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          do nn=1,ntot      

!----------------------------------------------------------------------
!    direct quantities                                            
!----------------------------------------------------------------------
            ww(1) = omegastr2(nn)
            ww(2) = gstr2(nn)
            ww(3) = taustr2(nn)
            ww(4) = cosangzk2(nn)

            qq(1)     = 3.0 * ( 1.0 - ww(1) )
            qq(2)         = 1.0 - ww(1) * ww(2)
            qq(3)     = qq(1)/qq(2)
            qq(4) = sqrt( qq(1) * qq(2) )
            qq(5) = sqrt (qq(3))
            qq(6) = 1.0 + twodi3 * qq(5)         
            qq(7) = 1.0 - twodi3 * qq(5)       

            rr(1) = 1./qq(6)
            rr(2) = qq(7)*rr(1)
            rr(3) = exp( -ww(3)          * qq(4) )
            rr(4) = 1.0/rr(3)
            rr(5) = 1.0/(qq(6) * rr(4) - qq(7) * rr(3) * rr(2) )

            ss(1) = 0.75 * ww(1)/(1.0 - (qq(4)*ww(4)      ) ** 2 )
            ss(2) = ss(1)*ww(4)*( 1.0 + ww(2)*qq(1)*onedi3)
            ss(3) = ss(1)*(1.0 + ww(2)*qq(1)*ww(4)** 2 )
            ss(4) = ss(2) - twodi3*ss(3)     
            ss(5) = ss(2) + twodi3 * ss(3)     
            ss(6) = exp( -ww(3)          / ww(4) )
            ss(7) = (ss(4)*ss(6) - ss(5)*rr(3)*rr(2))*rr(5)
            ss(8) = (ss(5) - qq(7)*ss(7))*rr(1)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
            rlayerdir2(nn) = qq(7) * ss(8) + qq(6)*ss(7) - ss(4)
            tlayerdir2(nn) = ((rr(3) * qq(6) * ss(8) + &
                               qq(7) * rr(4) * ss(7) -  &
                               ss(5) * ss(6) ) + ss(6) )
            tlayerde2(nn) = ss(6)

!----------------------------------------------------------------------c
!    diffuse quantities                                       
!    notes: the number of streams for the diffuse beam is fixed at 4.   
!    this calculation is done only for ng=1.                 
!----------------------------------------------------------------------c
            if (present (tlayerdif) .and. present(rlayerdif)) then
              if ( ng.eq.1 ) then   
                rsum = 0.0
                tsum = 0.0
                do ns = 1,NSTREAMS
                  tt(1) = 0.75 * ww(1)            / ( 1.0 - ( qq(4) * &
                          cosangstr(ns) ) ** 2 )
                  tt(2) = tt(1) * cosangstr(ns) * ( 1.0 +  &
                          ww(2)        * qq(1) * onedi3 )
                  tt(3) = tt(1) * ( 1.0 + ww(2)        * qq(1)*&
                          cosangstr(ns) ** 2 )
                  tt(4) = tt(2) - twodi3 * tt(3)
                  tt(5) = tt(2) + twodi3 * tt(3)
                  tt(6) = exp( -ww(3)          / cosangstr(ns) )
                  tt(7) = ( tt(4) * tt(6) - tt(5) *  &
                          rr(3) * rr(2)   ) * rr(5)
                  tt(8) = ( tt(5) - qq(7) * tt(7) )*rr(1)
                  if (nstr4) then
                    rsum = rsum + (qq(7)*tt(8) + qq(6)*tt(7) - tt(4))* &
                           wtstr(ns)*cosangstr(ns)
                    tsum = tsum + ((rr(3)*qq(6)*tt(8) +   &
                                    qq(7)*rr(4)*tt(7) -   &
                                    tt(5)*tt(6)) + tt(6))*  &
                                    wtstr(ns)*cosangstr(ns)
                  else 
                    rsum = rsum + (qq(7)*tt(8) + qq(6)*tt(7) - tt(4))
                    tsum = tsum + ( (rr(3)*qq(6)*tt(8) +    &
                                     qq(7)*rr(4)*tt(7) -   &
                                     tt(5)*tt(6)) + tt(6))
                  endif
                end do
                sumr(nn) = rsum
                sumt(nn) = tsum
              endif  !  ng == 1
            endif ! (present (tlayerdiff))
          end do  ! ntot loop

!---------------------------------------------------------------------
!     return results in proper locations in (i,j,k) arrays
!---------------------------------------------------------------------
          if (present(cloud)) then
            do i=1,ix           
              if (cloud(i,j,k) ) then
                rlayerdir(i,j,k) = rlayerdir2(cld_index(i,j,k))
                tlayerdir(i,j,k) = tlayerdir2(cld_index(i,j,k))
                tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
                if (present(tlayerdif)) then
                  if (ng .eq. 1) then
                    rlayerdif(i,j,k) = sumr(cld_index(i,j,k))
                    tlayerdif(i,j,k) = sumt(cld_index(i,j,k))
                  endif
                endif
              endif
            end do

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          else
            do i=1,ix            
              if (daylight(i,j) ) then
                rlayerdir(i,j,k) = rlayerdir2(cld_index(i,j,k))
                tlayerdir(i,j,k) = tlayerdir2(cld_index(i,j,k))
                tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
                if (present(tlayerdif)) then
                  if (ng .eq. 1) then
                    rlayerdif(i,j,k) = sumr(cld_index(i,j,k))
                    tlayerdif(i,j,k) = sumt(cld_index(i,j,k))
                  endif
                endif
              endif
            end do
          endif
        end do
      end do

!---------------------------------------------------------------------
 
end subroutine deledd

!#####################################################################
 subroutine deledd_4stream (ix, jx, kx, taustr, omegastr, pmstr1, pmstr2, pmstr3, cosang, ng, &
                   daylight, rlayerdir, tlayerdir, tlayerde,   &
                   rlayerdif, tlayerdif, cloud)
 implicit none
!---------------------------------------------------------------------- 
!    deledd calculates the reflection and transmission in the 
!    scattering layers using the delta-eddington method.         
!    references:                                                   
!      joseph, j.h., w. wiscombe, and j.a. weinman, the delta-eddington
!      approximation for radiative flux transfer.,j. atmos. sci.,33,  
!      2452-2459, 1976.                                              
!-------------------------------------------------------------------

integer,                   intent(in)              :: ix, jx, kx
real, dimension(:,:,:),    intent(inout)           :: taustr, omegastr
!stunew for 4 stream code
real, dimension(:,:,:),    intent(in)              :: pmstr1, pmstr2, pmstr3
!
real, dimension(:,:),    intent(in)                ::  cosang
integer,                   intent(in)              :: ng
logical, dimension(:,:),   intent(in)              :: daylight
real, dimension(:,:,:),    intent(out)             :: rlayerdir,   &
                                                      tlayerdir,   &
                                                      tlayerde
real, dimension(:,:,:),    intent(inout), optional :: rlayerdif,   &
                                                      tlayerdif
logical, dimension(:,:,:), intent(in), optional    :: cloud          

!----------------------------------------------------------------------
!  intent(in) variables:
!
!    ix,jx,kx
!    cosang      the cosine of the solar zenith angle    
!    ng          the number of gaussian angles to compute the diurnally 
!                averaged solar radiation (=1 unless lswg = true)       
!    daylight
!
!  intent(inout) variables:
!
!    taustr      the scaled extinction optical depth                    
!    omegastr    the scaled single-scattering albedo               
!
!  intent(out) variables:
!
!    rlayerdir   the layer reflectivity to a direct incident beam      
!    tlayerdir   the layer transmissivity to a direct incident beam   
!    rlayerdif   the layer reflectivity to a diffuse incident beam   
!    tlayerdif   the layer transmissivity to a diffuse incident beam
!    tlayerde    the layer transmissivity (non-scattered) to the direct 
!                incident beam                                       
!
! intent(in),optional:
!
!    cloud       flag for existence of a cloud (used only in 'ovc' mode)
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      real        :: rsum, tsum
      integer     :: k, i, ns, j, nn, ntot

      integer, dimension(ix, jx, kx) :: cld_index

      real,    dimension(ix)                  ::   &
                                          taustr2, omegastr2, &
                                           cosangzk2, rlayerdir2, &
                                           tlayerde2, tlayerdir2, &
                                           sumr, sumt
!stunew for 4 stream code
      real,    dimension(ix)                  ::   &
               pmstr12, pmstr22, pmstr32, rlayerdif2, tlayerdif2
      real        :: h(4),w(4),x(4),y(4),z(4),alambda(2)
      real        :: a(0:3),b(0:3),eta(0:3),tr(0:2)
      real        :: factor(2),p(2),q(2),r(2),c(2),d(2)
      real        :: rcossolar,delta,det

!----------------------------------------------------------------------
!  local variables:
!
!      i,j,k
!      ns
!      nn
!      ntot
!
!---------------------------------------------------------------------
!
!stunew define the gaussian weights and angles for the four stream diffuse
!calculation so that the setting of the number of streams is not considered
      real        :: WTSTR(4),cosangstr(4)

      wtstr(1)=0.347854845
      wtstr(2)=0.652145155
      wtstr(3)=0.347854845
      wtstr(4)=0.652145155
      cosangstr(1) = 0.069431844
      cosangstr(2) = 0.330009478
      cosangstr(3) = 0.930568156
      cosangstr(4) = 0.669990522
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do k=1,kx         
        do j=1,jx         

!---------------------------------------------------------------------
!    overcast sky mode. note: in this mode, the delta-eddington method
!    is performed only for spatial points containing a cloud.   
!-------------------------------------------------------------------
          nn = 0
          if (present(cloud)) then
            do i=1,ix          
              if (cloud(i,j,k) ) then
                nn = nn + 1
                cld_index(i,j,k) = nn
                taustr2(nn) = taustr(i,j,k)
                omegastr2(nn) = omegastr(i,j,k)
                cosangzk2(nn) = cosang(i,j)
!stunew for 4 stream code
                pmstr12(nn) = pmstr1(i,j,k)
                pmstr22(nn) = pmstr2(i,j,k)
                pmstr32(nn) = pmstr3(i,j,k)

!----------------------------------------------------------------------
!    note: the following are done to avoid the conservative scattering 
!    case, and to eliminate floating point errors in the exponential 
!    calculations, respectively.                      
!----------------------------------------------------------------------c
                if (omegastr2(nn) >= 9.9999999E-01) omegastr2(nn) = 9.9999999E-01
                if (taustr2(nn) >= 1.0E+02) taustr2(nn) = 1.0E+02
              endif
            end do

!----------------------------------------------------------------------c
!    clear sky mode. note: in this mode, the delta-eddington method is 
!    performed for all spatial points.                 
!----------------------------------------------------------------------c
          else
            do i=1,ix         
              if (daylight(i,j) ) then
                nn = nn + 1
                cld_index(i,j,k) = nn
                taustr2(nn) = taustr(i,j,k)
                omegastr2(nn) = omegastr(i,j,k)
                cosangzk2(nn) = cosang(i,j   )
!stunew for 4 stream code
                pmstr12(nn) = pmstr1(i,j,k)
                pmstr22(nn) = pmstr2(i,j,k)
                pmstr32(nn) = pmstr3(i,j,k)

!----------------------------------------------------------------------c
!    note: the following are done to avoid the conservative scattering  
!    case, and to eliminate floating point errors in the exponential 
!    calculations, respectively.                    
!----------------------------------------------------------------------c
                if (omegastr2(nn) >=9.9999999E-01) omegastr2(nn) = 9.9999999E-01
                if (taustr2(nn) >= 1.0E+02) taustr2(nn) = 1.0E+02
              endif
            end do
          endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          ntot = nn

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

          do nn=1,ntot      
! 
            a(0)=1.-omegastr2(nn)
            a(1)=3.-omegastr2(nn)*pmstr12(nn)
            a(2)=5.-omegastr2(nn)*pmstr22(nn)
            a(3)=7.-omegastr2(nn)*pmstr32(nn)
!  
            factor(1)=a(0)*a(1)+a(2)*a(3)/9.+4./9.*a(0)*a(3)
            factor(2)=1./9.*a(0)*a(1)*a(2)*a(3)
        !    write(3322,*) nn,ntot, a(0:4)
            alambda(1)=sqrt(0.5*(factor(1)+sqrt(factor(1)**2-4.*factor(2))))
            alambda(2)=sqrt(0.5*(factor(1)-sqrt(factor(1)**2-4.*factor(2))))
        !    write(3322,*) alambda(1),  alambda(2)
            tr(1)=exp(-alambda(1)*taustr2(nn))
            tr(2)=exp(-alambda(2)*taustr2(nn))
!
            p(1)=-a(0)/alambda(1)
            q(1)=0.5*(a(0)*a(1)/alambda(1)**2-1.)      
            r(1)=-1.5*(a(0)*a(1)/alambda(1)-alambda(1))/a(3) 
! 
            p(2)=-a(0)/alambda(2)
            q(2)=0.5*(a(0)*a(1)/alambda(2)**2-1.)     
            r(2)=-1.5*(a(0)*a(1)/alambda(2)-alambda(2))/a(3) 
!
            w(1)=(0.5-p(1)+0.625*q(1))
            w(2)=(0.5+p(1)+0.625*q(1))
            w(3)=(0.5-p(2)+0.625*q(2))
            w(4)=(0.5+p(2)+0.625*q(2))
!
            x(1)=(-0.125+0.625*q(1)-r(1))
            x(2)=(-0.125+0.625*q(1)+r(1))
            x(3)=(-0.125+0.625*q(2)-r(2))
            x(4)=(-0.125+0.625*q(2)+r(2))
!
            y(1)=(0.5+p(1)+0.625*q(1))*tr(1)
            y(2)=(0.5-p(1)+0.625*q(1))/tr(1)
            y(3)=(0.5+p(2)+0.625*q(2))*tr(2)
            y(4)=(0.5-p(2)+0.625*q(2))/tr(2)
!
            z(1)=(-0.125+0.625*q(1)+r(1))*tr(1)
            z(2)=(-0.125+0.625*q(1)-r(1))/tr(1)
            z(3)=(-0.125+0.625*q(2)+r(2))*tr(2)
            z(4)=(-0.125+0.625*q(2)-r(2))/tr(2)
!
            b(0)=omegastr2(nn)*0.25
            b(1)=-omegastr2(nn)*pmstr12(nn)*cosangzk2(nn)*0.25
            b(2)=omegastr2(nn)*pmstr22(nn)*(3.*cosangzk2(nn)**2-1.)*0.125
            b(3)=-omegastr2(nn)*pmstr32(nn)*(5.*cosangzk2(nn)**3-3.*cosangzk2(nn))*0.125
!  
            rcossolar=1./cosangzk2(nn) 
!
            tr(0)=exp(-taustr2(nn)*rcossolar)
!    
            delta=9.*rcossolar**4-rcossolar**2*(9.*a(0)*a(1)+a(2)*a(3)+4.*a(0)*a(3))+ &
                      a(0)*a(1)*a(2)*a(3)
            delta=1./delta
!
            eta(0)=((a(1)*b(0)-rcossolar*b(1))*(a(2)*a(3)-9.*rcossolar**2)+ &
              2.*rcossolar**2*(a(3)*b(2)-2.*a(3)*b(0)-3.*rcossolar*b(3)))*delta
!            
            eta(1)=((a(0)*b(1)-rcossolar*b(0))*(a(2)*a(3)-9.*rcossolar**2)- &
              2.*rcossolar*a(0)*(a(3)*b(2)-3.*rcossolar*b(3)))*delta
!
            eta(2)=((a(3)*b(2)-3.*rcossolar*b(3))*(a(0)*a(1)-rcossolar**2)- &
             2.*rcossolar*a(3)*(a(0)*b(1)-rcossolar*b(0)))*delta
!
            eta(3)=((a(2)*b(3)-3.*rcossolar*b(2))*(a(0)*a(1)-rcossolar**2)+ &
             rcossolar**2*(6.*a(0)*b(1)-4.*a(0)*b(3)-6.*rcossolar*b(0)))*delta
!
            h(1)=-(0.5*eta(0)-eta(1)+0.625*eta(2))
            h(2)=-(-0.125*eta(0)+0.625*eta(2)-eta(3))
            h(3)=-(0.5*eta(0)+eta(1)+0.625*eta(2))*tr(0)
            h(4)=-(-0.125*eta(0)+0.625*eta(2)+eta(3))*tr(0)
!
            det=(w(1)*x(2)-x(1)*w(2))*(y(3)*z(4)-z(3)*y(4))- &
              (w(1)*x(3)-x(1)*w(3))*(y(2)*z(4)-z(2)*y(4))+ &
              (w(1)*x(4)-x(1)*w(4))*(y(2)*z(3)-z(2)*y(3))+ &
              (w(2)*x(3)-x(2)*w(3))*(y(1)*z(4)-z(1)*y(4))- &
              (w(2)*x(4)-x(2)*w(4))*(y(1)*z(3)-z(1)*y(3))+ &
              (w(3)*x(4)-x(3)*w(4))*(y(1)*z(2)-z(1)*y(2)) 
            det=1./det
!
            c(1)=((h(1)*x(2)-h(2)*w(2))*(y(3)*z(4)-z(3)*y(4))- &
              (h(1)*x(3)-h(2)*w(3))*(y(2)*z(4)-z(2)*y(4))+ &
              (h(1)*x(4)-h(2)*w(4))*(y(2)*z(3)-z(2)*y(3))+ &
              (w(2)*x(3)-x(2)*w(3))*(h(3)*z(4)-h(4)*y(4))- &
              (w(2)*x(4)-x(2)*w(4))*(h(3)*z(3)-h(4)*y(3))+ &
              (w(3)*x(4)-x(3)*w(4))*(h(3)*z(2)-h(4)*y(2)))*det
!
            d(1)=((w(1)*h(2)-x(1)*h(1))*(y(3)*z(4)-z(3)*y(4))- &
              (w(1)*x(3)-x(1)*w(3))*(h(3)*z(4)-h(4)*y(4))+ &
              (w(1)*x(4)-x(1)*w(4))*(h(3)*z(3)-h(4)*y(3))+ &
              (h(1)*x(3)-h(2)*w(3))*(y(1)*z(4)-z(1)*y(4))- &
              (h(1)*x(4)-h(2)*w(4))*(y(1)*z(3)-z(1)*y(3))+ &
              (w(3)*x(4)-x(3)*w(4))*(y(1)*h(4)-z(1)*h(3)))*det 
!
            c(2)=((w(1)*x(2)-x(1)*w(2))*(h(3)*z(4)-h(4)*y(4))- &
              (w(1)*h(2)-x(1)*h(1))*(y(2)*z(4)-z(2)*y(4))+ &
              (w(1)*x(4)-x(1)*w(4))*(y(2)*h(4)-z(2)*h(3))+ &
              (w(2)*h(2)-x(2)*h(1))*(y(1)*z(4)-z(1)*y(4))- &
              (w(2)*x(4)-x(2)*w(4))*(y(1)*h(4)-z(1)*h(3))+ &
              (h(1)*x(4)-h(2)*w(4))*(y(1)*z(2)-z(1)*y(2)))*det
!
            d(2)=((w(1)*x(2)-x(1)*w(2))*(y(3)*h(4)-z(3)*h(3))- &
              (w(1)*x(3)-x(1)*w(3))*(y(2)*h(4)-z(2)*h(3))+ &
              (w(1)*h(2)-x(1)*h(1))*(y(2)*z(3)-z(2)*y(3))+ &
              (w(2)*x(3)-x(2)*w(3))*(y(1)*h(4)-z(1)*h(3))- &
              (w(2)*h(2)-x(2)*h(1))*(y(1)*z(3)-z(1)*y(3))+ &
              (w(3)*h(2)-x(3)*h(1))*(y(1)*z(2)-z(1)*y(2)))*det 
!
            rlayerdir2(nn)=((c(1)+d(1)+c(2)+d(2)+eta(0))+ &
              2.*(p(1)*(c(1)-d(1))+p(2)*(c(2)-d(2))+eta(1))+ &
              1.25*(q(1)*(c(1)+d(1))+q(2)*(c(2)+d(2))+ &
              eta(2)))*rcossolar
!
            tlayerdir2(nn)=((c(1)*tr(1)+d(1)/tr(1)+c(2)*tr(2)+ &
              d(2)/tr(2)+eta(0)*tr(0))-2.*(p(1)*(c(1)*tr(1)- &
              d(1)/tr(1))+p(2)*(c(2)*tr(2)-d(2)/tr(2))+ &
              eta(1)*tr(0))+1.25*(q(1)*(c(1)*tr(1)+d(1)/tr(1))+ &
              q(2)*(c(2)*tr(2)+d(2)/tr(2))+eta(2)*tr(0)))* &
              rcossolar+tr(0)
!
            tlayerde2(nn) = tr(0)
!
!----------------------------------------------------------------------c
!    diffuse quantities                                       
!    notes: the number of streams for the diffuse beam is fixed at 4.   
!    this calculation is done only for ng=4.                 
!----------------------------------------------------------------------c
            if (present (tlayerdif) .and. present(rlayerdif)) then

              rsum = 0.0
              tsum = 0.0

              do ns = 1,4
!
                b(1)=-omegastr2(nn)*pmstr12(nn)*cosangstr(ns)*0.25
                b(2)=omegastr2(nn)*pmstr22(nn)*(3.*cosangstr(ns)**2-1.)*0.125
                b(3)=-omegastr2(nn)*pmstr32(nn)*(5.*cosangstr(ns)**3-3.*cosangstr(ns))*0.125
!
                rcossolar=1./cosangstr(ns)
!
                tr(0)=exp(-taustr2(nn)*rcossolar)
!
                delta=9.*rcossolar**4-rcossolar**2*(9.*a(0)*a(1)+a(2)*a(3)+4.*a(0)*a(3))+ &
                      a(0)*a(1)*a(2)*a(3)
                delta=1./delta
!
                eta(0)=((a(1)*b(0)-rcossolar*b(1))*(a(2)*a(3)-9.*rcossolar**2)+ &
                  2.*rcossolar**2*(a(3)*b(2)-2.*a(3)*b(0)-3.*rcossolar*b(3)))*delta
!
                eta(1)=((a(0)*b(1)-rcossolar*b(0))*(a(2)*a(3)-9.*rcossolar**2)- &
                  2.*rcossolar*a(0)*(a(3)*b(2)-3.*rcossolar*b(3)))*delta
!
                eta(2)=((a(3)*b(2)-3.*rcossolar*b(3))*(a(0)*a(1)-rcossolar**2)- &
                  2.*rcossolar*a(3)*(a(0)*b(1)-rcossolar*b(0)))*delta
!
                eta(3)=((a(2)*b(3)-3.*rcossolar*b(2))*(a(0)*a(1)-rcossolar**2)+ &
                    rcossolar**2*(6.*a(0)*b(1)-4.*a(0)*b(3)-6.*rcossolar*b(0)))*delta
!
                h(1)=-(0.5*eta(0)-eta(1)+0.625*eta(2))
                h(2)=-(-0.125*eta(0)+0.625*eta(2)-eta(3))
                h(3)=-(0.5*eta(0)+eta(1)+0.625*eta(2))*tr(0)
                h(4)=-(-0.125*eta(0)+0.625*eta(2)+eta(3))*tr(0)
!
                det=(w(1)*x(2)-x(1)*w(2))*(y(3)*z(4)-z(3)*y(4))- &
                  (w(1)*x(3)-x(1)*w(3))*(y(2)*z(4)-z(2)*y(4))+ &
                  (w(1)*x(4)-x(1)*w(4))*(y(2)*z(3)-z(2)*y(3))+ &
                  (w(2)*x(3)-x(2)*w(3))*(y(1)*z(4)-z(1)*y(4))- &
                  (w(2)*x(4)-x(2)*w(4))*(y(1)*z(3)-z(1)*y(3))+ &
                  (w(3)*x(4)-x(3)*w(4))*(y(1)*z(2)-z(1)*y(2))
                det=1./det

                c(1)=((h(1)*x(2)-h(2)*w(2))*(y(3)*z(4)-z(3)*y(4))- &
                  (h(1)*x(3)-h(2)*w(3))*(y(2)*z(4)-z(2)*y(4))+ &
                  (h(1)*x(4)-h(2)*w(4))*(y(2)*z(3)-z(2)*y(3))+ &
                  (w(2)*x(3)-x(2)*w(3))*(h(3)*z(4)-h(4)*y(4))- &
                  (w(2)*x(4)-x(2)*w(4))*(h(3)*z(3)-h(4)*y(3))+ &
                  (w(3)*x(4)-x(3)*w(4))*(h(3)*z(2)-h(4)*y(2)))*det
!
                d(1)=((w(1)*h(2)-x(1)*h(1))*(y(3)*z(4)-z(3)*y(4))- &
                  (w(1)*x(3)-x(1)*w(3))*(h(3)*z(4)-h(4)*y(4))+ &
                  (w(1)*x(4)-x(1)*w(4))*(h(3)*z(3)-h(4)*y(3))+ &
                  (h(1)*x(3)-h(2)*w(3))*(y(1)*z(4)-z(1)*y(4))- &
                  (h(1)*x(4)-h(2)*w(4))*(y(1)*z(3)-z(1)*y(3))+ &
                  (w(3)*x(4)-x(3)*w(4))*(y(1)*h(4)-z(1)*h(3)))*det
!
                c(2)=((w(1)*x(2)-x(1)*w(2))*(h(3)*z(4)-h(4)*y(4))- &
                  (w(1)*h(2)-x(1)*h(1))*(y(2)*z(4)-z(2)*y(4))+ &
                  (w(1)*x(4)-x(1)*w(4))*(y(2)*h(4)-z(2)*h(3))+ &
                  (w(2)*h(2)-x(2)*h(1))*(y(1)*z(4)-z(1)*y(4))- &
                  (w(2)*x(4)-x(2)*w(4))*(y(1)*h(4)-z(1)*h(3))+ &
                  (h(1)*x(4)-h(2)*w(4))*(y(1)*z(2)-z(1)*y(2)))*det
!
                d(2)=((w(1)*x(2)-x(1)*w(2))*(y(3)*h(4)-z(3)*h(3))- &
                  (w(1)*x(3)-x(1)*w(3))*(y(2)*h(4)-z(2)*h(3))+ &
                  (w(1)*h(2)-x(1)*h(1))*(y(2)*z(3)-z(2)*y(3))+ &
                  (w(2)*x(3)-x(2)*w(3))*(y(1)*h(4)-z(1)*h(3))- &
                  (w(2)*h(2)-x(2)*h(1))*(y(1)*z(3)-z(1)*y(3))+ &
                  (w(3)*h(2)-x(3)*h(1))*(y(1)*z(2)-z(1)*y(2)))*det

                  rsum = rsum + (((c(1)+d(1)+c(2)+d(2)+eta(0))+ &
                                  2.*(p(1)*(c(1)-d(1))+p(2)*(c(2)-d(2))+eta(1))+ &
                                  1.25*(q(1)*(c(1)+d(1))+q(2)*(c(2)+d(2))+eta(2)))*rcossolar)*cosangstr(ns)*wtstr(ns)
                  tsum = tsum + (((c(1)*tr(1)+d(1)/tr(1)+c(2)*tr(2)+d(2)/tr(2)+eta(0)*tr(0))- &
                                  2.*(p(1)*(c(1)*tr(1)-d(1)/tr(1))+p(2)*(c(2)*tr(2)- &
                                  d(2)/tr(2))+eta(1)*tr(0))+1.25*(q(1)*(c(1)*tr(1)+d(1)/tr(1))+ &
                                  q(2)*(c(2)*tr(2)+d(2)/tr(2))+eta(2)*tr(0)))*rcossolar+tr(0))*cosangstr(ns)*wtstr(ns)
            end do

            sumr(nn) = rsum
            sumt(nn) = tsum

            endif ! (present (tlayerdiff))
!
          end do  ! ntot loop

!---------------------------------------------------------------------
!     return results in proper locations in (i,j,k) arrays
!---------------------------------------------------------------------
          if (present(cloud)) then
            do i=1,ix           
              if (cloud(i,j,k) ) then
                rlayerdir(i,j,k) =rlayerdir2(cld_index(i,j,k))
                tlayerdir(i,j,k) =tlayerdir2(cld_index(i,j,k))
                tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
                rlayerdif(i,j,k) =sumr(cld_index(i,j,k))
                tlayerdif(i,j,k) =sumt(cld_index(i,j,k))
              endif
            end do

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          else
            do i=1,ix            
              if (daylight(i,j) ) then
                rlayerdir(i,j,k) =rlayerdir2(cld_index(i,j,k))
                tlayerdir(i,j,k) =tlayerdir2(cld_index(i,j,k))
                tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
                rlayerdif(i,j,k) =sumr(cld_index(i,j,k))
                tlayerdif(i,j,k) =sumt(cld_index(i,j,k))
              endif
            end do
          endif
        end do
      end do

!---------------------------------------------------------------------
 
end subroutine deledd_4stream

                   end module esfsw_driver_mod

