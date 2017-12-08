                 module esfsw_parameters_mod
!
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Code to initialize shortwave parameters and access flux data.
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes shortwave radiation calculation parameters such as
!  solar flux input at top of the atmosphere, number of shortwave bands
!  depending on the spectral resolution used, number of frequency points
!  in the gaussian quadrature algorithm, the number of streams used in
!  multiple stream flux algorithm, and the number of water vapor bands.
!
!  The code also provides two access methods: get and put solar flux data
!  
! </DESCRIPTION>

!    shared modules:

use mpp_mod,           only: input_nml_file
use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, close_file

!--------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!     esfsw_parameters_mod defines parameters for esf shortwave code,
!     including a description of the band structure  used to define the
!     solar spectrum.
!------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'

!--------------------------------------------------------------------
!----- interfaces ------

public       &
         esfsw_parameters_init, &
         esfsw_parameters_end

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)  :: sw_resolution = 'low' ! no longer does anything
                                            ! superceeded by sw_code_version
character(len=16)  :: sw_code_version = 'esf1999' ! define version of code
                                            ! either esf1999 or esf2015
integer            :: sw_bands = 18         ! Number of shortwave
                                            ! bands only works if 
                                            ! sw_code_version = esf2015
integer            :: sw_int_points = 74    ! number of shortwave
                                            ! frequency intergration points only works
                                            ! if sw_code_version = esf2015
integer            :: sw_NIRVISgas_bands = 11   ! number of shortwave
                                            ! bands in the NIR and VIS only works
                                            ! if sw_code_version = esf2015
integer            :: sw_diff_streams = 0   ! number of streams of
                                            ! diffuse radiation that
                                            ! are considered


namelist /esfsw_parameters_nml/    &
                                 sw_resolution,   &
                                 sw_diff_streams, &
                                 sw_code_version, &
                                 sw_bands,        &
                                 sw_int_points,   &
                                 sw_NIRVISgas_bands

!-------------------------------------------------------------------
!----- public data --------

!---------------------------------------------------------------------
!    TOT_WVNUMS_LOCAL  number of wavenumbers included in the
!                      parameterization of the solar spectrum 
!DEL Solar_spect    solar_spectrum_type variable defining the nature
!DEL                of the solar spectral paramaterization
!---------------------------------------------------------------------
integer, parameter :: TOT_WVNUMS_LOCAL  = 57600


!-------------------------------------------------------------------
!----- private data --------

logical :: module_is_initialized = .false.  ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                      contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     
!                     PUBLIC SUBROUTINES
!            
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
!
! <SUBROUTINE NAME="esfsw_parameters_init">
!  <OVERVIEW>
!   Subroutine that initializes and set up shortwave radiation.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine that initializes shortwave radiation calculation parameters such as
!   solar flux input at top of the atmosphere, number of shortwave bands
!   depending on the spectral resolution used, number of frequency points
!   in the gaussian quadrature algorithm, the number of streams used in
!   multiple stream flux algorithm, and the number of water vapor bands.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call esfsw_parameters_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine esfsw_parameters_init (nbands, nfrqpts, nh2obands, &
                                  nstreams, tot_wvnums)

integer, intent(out),target :: nbands, nfrqpts, nh2obands, &
                               nstreams, tot_wvnums

!------------------------------------------------------------------
!    esfsw_parameters_init is the constructor for esfsw_parameters_mod.
!------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      integer    ::  unit, ierr, io, logunit

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
      read (input_nml_file, nml=esfsw_parameters_nml, iostat=io)
      ierr = check_nml_error(io,'esfsw_parameters_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=esfsw_parameters_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'esfsw_parameters_nml')
        end do
10      call close_file (unit)
      endif
#endif

!--------------------------------------------------------------------
!    process the namelist entries to obtain the parameters specifying
!    the solar spectral parameterization.
!--------------------------------------------------------------------
      if (trim(sw_code_version) == 'esf1999') then
        nbands = 18
        nfrqpts = 38
        nh2obands = 9
      else if (trim(sw_code_version) == 'esf2015') then
        nbands = sw_bands 
        nfrqpts = sw_int_points 
        nh2obands = sw_NIRVISgas_bands      
      else       
        call error_mesg ( 'esfsw_parameters_mod',   &
           'esf_version must be specified as "esf2015" or "es1999', FATAL)
      endif

      if (sw_diff_streams == 4) then
        nstreams = 4
      else if (sw_diff_streams == 1) then
        nstreams = 1
      else
        call error_mesg ( 'esfsw_parameters_mod',   &
          ' sw_diff_streams must be specified as either 1 or 4.', FATAL)
      endif

!---------------------------------------------------------------------
!    include the total number of wavenumbers in the solar parameter-
!    ization in the solar_spectrum_type variable.
!---------------------------------------------------------------------
      tot_wvnums = TOT_WVNUMS_LOCAL

!---------------------------------------------------------------------
!    write version number and namelist to logfile also write out
!    some key parameters obtained from an input data file.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) then
        write (logunit,9000) nbands, nfrqpts, nstreams, nh2obands, tot_wvnums
        write (logunit, nml=esfsw_parameters_nml)
      endif  

!------------------------------------------------------------------
!    mark the module as initialized.
!------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------
9000  format ( '  NBANDS=  ', i4, '  NFRQPTS=', i4, &
               '  NSTREAMS= ', i4, '  NH2OBANDS= ', i4, &
               '  TOT_WVNUMS= ',i5 )

!------------------------------------------------------------------


end subroutine esfsw_parameters_init



!####################################################################
!
! <SUBROUTINE NAME="esfsw_parameters_end">
!  <OVERVIEW>
!   Subroutine that is the destructor for esfsw_parameters_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine that deallocates module variables and marks the module  
!   as uninitialized.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call esfsw_parameters_end
!  </TEMPLATE>
! </SUBROUTINE>
!

subroutine esfsw_parameters_end

!--------------------------------------------------------------------
!    esfsw_parameters_end is the destructor for esfsw_parameters_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('esfsw_parameters_mod',   &
             'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------

end subroutine esfsw_parameters_end

!###################################################################

      end module esfsw_parameters_mod

