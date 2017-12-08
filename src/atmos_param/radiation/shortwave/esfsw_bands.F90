
module esfsw_bands_mod

use fms_mod,      only: write_version_number, &
                        error_mesg, FATAL

use esfsw_utilities_mod, only: thickavg, thinavg

!--------------------------------------------------------------------

implicit none 
private

!---------------------------------------------------------------------
!-------  version number --------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------

public :: esfsw_bands_init, &
          esfsw_bands_end, &
          esfsw_number_of_bands, &
          esfsw_bands, &
          esfsw_band_segments, &
          esfsw_solar_flux_init, &
          esfsw_solar_flux, &
          esfsw_thickavg, &
          esfsw_thinavg

interface esfsw_thickavg
   module procedure esfsw_thickavg_3d
   module procedure esfsw_thickavg_0d
   module procedure esfsw_thickavg_1band
   module procedure esfsw_thickavg_isccp
end interface

!---------------------------------------------------------------------
!------- private data ------

! wave numbers in 1/cm
real    :: wvnum_340 = 1.0E+04/0.34        ! near uv
real    :: wvnum_380 = 1.0E+04/0.38
real    :: wvnum_440 = 1.0E+04/0.44        ! blue
real    :: wvnum_550 = 1.0E+04/0.55        ! visible light (0.55 microns)
real    :: wvnum_670 = 1.0E+04/0.67        ! red
real    :: wvnum_870 = 1.0E+04/0.87
real    :: wvnum_one_micron = 1.0E+04/1.00 ! near ir
real    :: wvnum_onepsix_micron = 1.0E+04/1.61

!---------------------------------------------------------------------
!--- module variables ---
!
!    endwvnbands       Highest wave number in each of the solar
!                      spectral parameterization bands [cm(-1)].
!                      This array has a starting index of zero.
!
!    solflxbandref     Solar flux in each parameterization band,
!                      used for defining band-averaged optical parameters.
!                      If the solar constant is time-invariant,
!                      it is also the solar flux in each
!                      parameterization band (solflxband).
!                                
!    solarfluxtoa      Highly-resolved solar flux at toa in 
!                      "tot_wvnums" bands
!
!---------------------------------------------------------------------

integer, allocatable :: endwvnbands_local(:)
real,    allocatable :: solarfluxtoa_local(:)
real,    allocatable :: solflxbandref_local(:)
real,    allocatable :: solflxband_local(:)
real    :: solar_constant_local
logical :: solflxband_initialized = .false.
integer :: nbands, tot_wvnums

logical :: module_is_initialized=.false.   ! module is initialized ?


CONTAINS

!#####################################################################
!
!                     PUBLIC SUBROUTINES
!
!#####################################################################

  subroutine esfsw_bands_init ( endwvnbands, solflxbandref, solarfluxtoa )
  integer, intent(in) :: endwvnbands(0:)
  real,    intent(in) :: solflxbandref(:), solarfluxtoa(:)

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    write version number to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)

!---------------------------------------------------------------------
!    save the shortwave band information
!---------------------------------------------------------------------
      nbands = size(endwvnbands,1)-1
      tot_wvnums = size(solarfluxtoa,1)

      allocate(endwvnbands_local(0:nbands))
      allocate(solflxbandref_local(nbands))
      allocate(solflxband_local(nbands))
      allocate(solarfluxtoa_local(tot_wvnums))
      endwvnbands_local = endwvnbands
      solflxbandref_local = solflxbandref
      solarfluxtoa_local = solarfluxtoa

!-------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------

  end subroutine esfsw_bands_init

!#####################################################################

 subroutine esfsw_number_of_bands (nbands)
 integer, intent(out) :: nbands

!---------------------------------------------------------------------
!    make sure module has been initialized
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then 
        call error_mesg ('esfsw_bands_mod',   &
              'module has not been initialized', FATAL )
      endif
!---------------------------------------------------------------------

      nbands = size(solflxbandref_local,1)
      
!---------------------------------------------------------------------

 end subroutine esfsw_number_of_bands

!#####################################################################

subroutine esfsw_solar_flux_init (solcon, solflxband) 
real, intent(in) :: solcon, solflxband(:)

!---------------------------------------------------------------------
!    make sure module has been initialized
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then 
        call error_mesg ('esfsw_bands_mod',   &
              'module has not been initialized', FATAL )
      endif
!---------------------------------------------------------------------
      if (size(solflxband,1) .ne. size(solflxband_local,1)) then 
        call error_mesg ('esfsw_bands_mod',   &
             'incorrect number of bands in esfsw_solar_flux', FATAL)
      endif
!---------------------------------------------------------------------

      solar_constant_local = solcon
      solflxband_local = solflxband
      solflxband_initialized = .true.

!---------------------------------------------------------------------

end subroutine esfsw_solar_flux_init

!#####################################################################
!
! Returns the solar flux by band
!

subroutine esfsw_solar_flux (solflx, ref) 
real,    intent(out)          :: solflx(:)
logical, intent(in), optional :: ref
!---------------------------------------------------------------------
!    make sure module has been initialized
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then 
        call error_mesg ('esfsw_bands_mod',   &
              'module has not been initialized', FATAL )
      endif
!---------------------------------------------------------------------
      if (size(solflx,1) .ne. size(solflxbandref_local,1)) then 
        call error_mesg ('esfsw_bands_mod',   &
             'incorrect number of bands in esfsw_solar_flux', FATAL)
      endif
!---------------------------------------------------------------------
      if (present(ref)) then 
         if (ref) then 
            solflx = solflxbandref_local
            return
         endif
      endif

      if (solflxband_initialized) then 
         solflx = solflxband_local
      else 
         call error_mesg ('esfsw_bands_mod',   &
             'time varying solar flux not initialized', FATAL )
      endif
!---------------------------------------------------------------------

end subroutine esfsw_solar_flux

!#####################################################################

  subroutine esfsw_bands ( w340_band_indx, w380_band_indx, &
                           w440_band_indx, w550_band_indx, &
                           w670_band_indx, w870_band_indx, &
                           one_micron_indx, onepsix_micron_indx )

  integer, optional,      intent(out) :: w340_band_indx, &
                                         w380_band_indx, &
                                         w440_band_indx, &
                                         w550_band_indx, &
                                         w670_band_indx, &
                                         w870_band_indx, &
                                         one_micron_indx, &
                                         onepsix_micron_indx

  integer :: ni

!---------------------------------------------------------------------
!    make sure module has been initialized
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then 
        call error_mesg ('esfsw_bands_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    define the band index corresponding to near ultraviolet light
!    (0.34 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      if (present(w340_band_indx)) then
        do ni=1,nbands
          if (endwvnbands_local(ni) > wvnum_340) then
            w340_band_indx = ni
            exit
          endif
        end do
      endif
!---------------------------------------------------------------------
!    define the band index corresponding to near ultraviolet light
!    (0.38 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      if (present(w380_band_indx)) then
        do ni=1,nbands
          if (endwvnbands_local(ni) > wvnum_380) then
            w380_band_indx = ni
            exit
          endif
        end do
      endif
!---------------------------------------------------------------------
!    define the band index corresponding to blue light
!    (0.44 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      if (present(w440_band_indx)) then
        do ni=1,nbands
          if (endwvnbands_local(ni) > wvnum_440) then
            w440_band_indx = ni
            exit
          endif
        end do
      endif
!---------------------------------------------------------------------
!    define the band index corresponding to visible light
!    (0.55 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      if (present(w550_band_indx)) then
        do ni=1,nbands
          if (endwvnbands_local(ni) > wvnum_550) then
            w550_band_indx = ni
            exit
          endif
        end do
      endif
!---------------------------------------------------------------------
!    define the band index corresponding to red light
!    (0.67 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      if (present(w670_band_indx)) then
        do ni=1,nbands
          if (endwvnbands_local(ni) > wvnum_670) then
            w670_band_indx = ni
            exit
          endif
        end do
      endif
!---------------------------------------------------------------------
!    define the band index corresponding to 870nm
!    (0.87 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      if (present(w870_band_indx)) then
        do ni=1,nbands
          if (endwvnbands_local(ni) > wvnum_870) then
            w870_band_indx = ni
            exit
          endif
        end do
      endif
!---------------------------------------------------------------------
!    define the band index corresponding to near infra red band
!    (1.00 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      if (present(one_micron_indx)) then
        do ni=1,nbands
          if (endwvnbands_local(ni) > wvnum_one_micron) then
            one_micron_indx = ni
            exit
          endif
        end do
      endif
!---------------------------------------------------------------------
!    define the band index corresponding to 
!    (1.61 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      if (present(onepsix_micron_indx)) then
        do ni=1,nbands
          if (endwvnbands_local(ni) > wvnum_onepsix_micron) then
            onepsix_micron_indx = ni
            exit
          endif
        end do
      endif


  end subroutine esfsw_bands

!#####################################################################

  subroutine esfsw_band_segments ( endwvn, solflx, ind1, ind2 )
  integer, intent(in)  :: endwvn(:)
  real,    intent(out) :: solflx(:,:)
  integer, intent(out) :: ind1(:), ind2(:)

  real :: sumsol
  integer :: ni, nb, nw

!---------------------------------------------------------------------
!    make sure module has been initialized
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then 
        call error_mesg ('esfsw_bands_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    make sure argument array sizes match initialized array sizes
!---------------------------------------------------------------------

      if (size(solflx,1) /= nbands .or. size(solflx,2) /= size(endwvn,1)) &
          call error_mesg( 'esfsw_bands_mod', &
                           'incorrect size for input arguments '//&
                           'to esfsw_band_segments', FATAL)

!--------------------------------------------------------------------
!    verify consistency between the highest wavenumber in the solar 
!    spectrum and the highest wavenumber in the various particle
!    spectral intervals to assure that all solar spectral bands are 
!    assigned to a particle band.
!--------------------------------------------------------------------

      if (tot_wvnums /= endwvn(size(endwvn,1))) &
          call error_mesg( 'esfsw_bands_mod', &
               'in subroutine esfsw_band_segments the '//&
               'highest wave number in particle spectrum '//&
               'differs from highest wavenumber in solar spectrum', FATAL)

!---------------------------------------------------------------------
!    the single scattering properties for wavenumbers < the lower limit
!    to which the various parameterizations apply are assigned the
!    values in the lowest interval of the parameterization; thus all
!    solar spectral parameterization bands are assigned to a 
!    water-species parameterization band.
!---------------------------------------------------------------------

      solflx(:,:) = 0.0
      ind1(1) = 1

      ni = 1
      nb = 1
      sumsol = 0.0
!---------------------------------------------------------------------
!    integrate over wavenumber, summing the solar spectral bands con-
!    sistent with the parameterization band structure and the particle 
!    band structure. when the loop is ended, the solar flux in each
!    spectral region generated by overlapping the parameterization
!    spectrum and a particular particle spectrum will be resident in
!    the solflx arrays.
!---------------------------------------------------------------------
      do nw = 1, tot_wvnums
         sumsol = sumsol + solarfluxtoa_local(nw)

         if ( nw == endwvn(ni) ) then
            solflx(nb,ni) = sumsol
            sumsol = 0.0
         end if

         if ( nw == endwvnbands_local(nb) ) then
            if ( nw /= endwvn(ni) ) then
               solflx(nb,ni) = sumsol
               sumsol = 0.0
            end if

            ind2(nb) = ni
            nb = nb+1
            if ( nb <= nbands ) then
               if ( nw == endwvn(ni) ) then
                  ind1(nb) = ni+1
               else
                  ind1(nb) = ni
               end if
            end if
         end if

         if ( nw == endwvn(ni) ) ni = ni+1
      end do

  end subroutine esfsw_band_segments

!#####################################################################
!
! <SUBROUTINE NAME="esfsw_thickavg_3d">
!  <OVERVIEW>
!   Subroutine to use thick-averaging technique to define band interval
!   single scattering properties.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thick-averaging technique to define the single-scattering    
! properties of the parameterization band spectral intervals from the  
! specified spectral intervals of the particular scatterer.            
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine esfsw_thickavg_3d (nivl1    , nivl2     , nivls   ,   &
!                        nbands, $
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        extband  , ssalbband , asymmband)
!  </TEMPLATE>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="nivls" TYPE="integer">
!   number of specified scattering spectral intervals
!  </IN>
!  <IN NAME="nbands" TYPE="integer">
!   number of spectral bands
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <INOUT NAME="ssalbivl" TYPE="real">
!   the specified spectral values of the single-scattering albedo
!  </INOUT>
!  <IN NAME="asymmivl" TYPE="real">
!   the specified spectral values of the asymmetry factor
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
!  <OUT NAME="ssalbband" TYPE="real">
!   the parameterization band values of the single-scattering albedo
!  </OUT>
!  <OUT NAME="asymmband" TYPE="real">
!   the parameterization band values of the asymmetry factor
!  </OUT>
! </SUBROUTINE>
!
subroutine esfsw_thickavg_3d (nivl1, nivl2, nivls, nbands, extivl, ssalbivl,&
                        asymmivl, solflxivl, mask, extband,&
                        ssalbband, asymmband)

!---------------------------------------------------------------------
!    esfsw_thickavg_3d uses the thick-averaging technique to define the 
!    single-scattering properties of the parameterization band spectral
!    intervals from the  specified spectral intervals of the particular
!    scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl
real, dimension(:,:,:,:), intent(out)      :: extband, ssalbband,   &
                                              asymmband
logical, dimension(:,:,:), intent(in)      :: mask

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    nivls       number of specified scattering spectral intervals      
!    nbands
!    extivl      specified spectral values of the extinction coefficient
!    asymmivl    the specified spectral values of the asymmetry     
!                factor                                           
!    solflxivl   the solar flux in each specified scattering spectral
!                interval                                         
!
!  intent(inout) variables:
!
!    ssalbivl    the specified spectral values of the single-       
!                scattering albedo                                   
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    ssalbband   the parameterization band values of the single-   
!                scattering albedo                                  
!    asymmband   the parameterization band values of the asymmetry   
!                factor                                               
!    
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!   module variable:
!
!    solflxbandref_local  the solar flux in each parameterization band  
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_bands_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

      if (nbands /= size(solflxbandref_local,1)) then
        call error_mesg ('esfsw_bands_mod',  &
         'nbands argument not equal to number of SW bands registered', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------

     call thickavg (nivl1, nivl2, nivls, nbands, extivl, ssalbivl,&
                    asymmivl, solflxivl, solflxbandref_local, &
                    mask, extband, ssalbband, asymmband)

!---------------------------------------------------------------------

end subroutine esfsw_thickavg_3d

!####################################################################
! <SUBROUTINE NAME="esfsw_thickavg_0d">
!  <OVERVIEW>
!   Subroutine to use thick-averaging technique to define band interval
!   single scattering properties.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thick-averaging technique to define the single-scattering    
! properties of the parameterization band spectral intervals from the  
! specified spectral intervals of the particular scatterer.            
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine esfsw_thickavg_0d (nivl1    , nivl2     , nivls   ,   &
!                        nbands,  &
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        extband  , ssalbband , asymmband)
!  </TEMPLATE>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="nivls" TYPE="integer">
!   number of specified scattering spectral intervals
!  </IN>
!  <IN NAME="nbands" TYPE="integer">
!   number of spectral bands
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <INOUT NAME="ssalbivl" TYPE="real">
!   the specified spectral values of the single-scattering albedo
!  </INOUT>
!  <IN NAME="asymmivl" TYPE="real">
!   the specified spectral values of the asymmetry factor
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
!  <OUT NAME="ssalbband" TYPE="real">
!   the parameterization band values of the single-scattering albedo
!  </OUT>
!  <OUT NAME="asymmband" TYPE="real">
!   the parameterization band values of the asymmetry factor
!  </OUT>
! </SUBROUTINE>
!
subroutine esfsw_thickavg_0d (nivl1, nivl2, nivls, nbands, extivl, &
                              ssalbivl, asymmivl, solflxivl, &
                              extband, ssalbband , asymmband)

!---------------------------------------------------------------------
!    esfsw_thickavg_0d uses the thick-averaging technique to define the 
!    single-scattering properties of the parameterization band spectral
!    intervals from the  specified spectral intervals of the particular
!    scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:),       intent(in)       :: extivl, asymmivl
real, dimension(:),       intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl
real, dimension(:),       intent(out)      :: extband, ssalbband, &
                                              asymmband

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    nivls       number of specified scattering spectral intervals      
!    nbands
!    extivl      specified spectral values of the extinction coefficient
!    asymmivl    the specified spectral values of the asymmetry     
!                factor                                           
!    solflxivl   the solar flux in each specified scattering spectral
!                interval                                         
!
!  intent(inout) variables:
!
!    ssalbivl    the specified spectral values of the single-       
!                scattering albedo                                   
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    ssalbband   the parameterization band values of the single-   
!                scattering albedo                                  
!    asymmband   the parameterization band values of the asymmetry   
!                factor                                               
!    
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  module variables:
!
!    solflxbandref_local  the solar flux in each parameterization band  
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_bands_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

      if (nbands /= size(solflxbandref_local,1)) then
        call error_mesg ('esfsw_bands_mod',  &
         'nbands argument not equal to number of SW bands registered', &
                                                                 FATAL)
      endif

!----------------------------------------------------------------------
! call thickavg_0d with solflxbandref

      call thickavg (nivl1, nivl2, nivls, nbands, extivl, ssalbivl, &
                     asymmivl, solflxivl, solflxbandref_local, &
                     extband, ssalbband , asymmband)
!----------------------------------------------------------------------

end subroutine esfsw_thickavg_0d

!##################################################################
!
! <SUBROUTINE NAME="esfsw_thickavg_isccp">
!  <OVERVIEW>
!   Subroutine to use thick-averaging technique to define band interval
!   single scattering properties for a single specified band.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thick-averaging technique to define the single-scattering    
! properties of the specified parameterization band spectral interval 
! from the specified spectral intervals of the particular scatterer.    
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine esfsw_thickavg (nband, nivl1, nivl2, extivl, solflxivl, &
!                             mask, extband )
!  </TEMPLATE>
!  <IN NAME="nband" TYPE="integer">
!
!  </IN>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <IN NAME="mask" TYPE="logical">
!   mask is .true. at gridpoints where extband needs to be calculated
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
! </SUBROUTINE>
!
subroutine esfsw_thickavg_isccp (nband, nivl1, nivl2, extivl,  &
                                 solflxivl, mask, extband)

!---------------------------------------------------------------------
!    esfsw_thickavg_isccp uses the thick-averaging technique to define the 
!    solar extinction for the single specified parameterization band 
!    spectral interval (nband) from the  specified spectral intervals 
!    of the particular scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer,                  intent(in)       :: nband
integer,                  intent(in)       :: nivl1, nivl2
real, dimension(:,:,:,:), intent(in)       :: extivl
real, dimension(:,:),     intent(in)       :: solflxivl
logical, dimension(:,:,:),intent(in)       :: mask
real, dimension(:,:,:),   intent(out)      :: extband

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nband       the sw parameterization band for which the optical
!                properties are being calculated
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    extivl      specified spectral values of the extinction coefficient
!    mask        logical indicating the points at which the band values 
!                should be calculated       
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    
!--------------------------------------------------------------------
!  module variables:
!
!    solflxbandref_local  the solar flux in each parameterization band  
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_bands_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
! call thickavg_isccp with solar band

      call thickavg (nband, nivl1, nivl2, extivl, solflxivl, &
                     solflxbandref_local(nband), mask, extband)
!--------------------------------------------------------------------

end subroutine esfsw_thickavg_isccp

!##################################################################
!
! <SUBROUTINE NAME="esfsw_thickavg_1band">
!  <OVERVIEW>
!   Subroutine to use thick-averaging technique to define band interval
!   single scattering properties for a single specified band.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thick-averaging technique to define the single-scattering    
! properties of the specified parameterization band spectral interval 
! from the specified spectral intervals of the particular scatterer.    
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine esfsw_thickavg (nband, nivl1  , nivl2   , nivls   , &
!                        nbands, &
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        mask, extband  , ssalbband ,  &
!                        asymmband)
!  </TEMPLATE>
!  <IN NAME="nband" TYPE="integer">
!
!  </IN>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="nivls" TYPE="integer">
!   number of specified scattering spectral intervals
!  </IN>
!  <IN NAME="nbands" TYPE="integer">
!   number of spectral bands
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <INOUT NAME="ssalbivl" TYPE="real">
!   the specified spectral values of the single-scattering albedo
!  </INOUT>
!  <IN NAME="asymmivl" TYPE="real">
!   the specified spectral values of the asymmetry factor
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <IN NAME="mask" TYPE="logical">
!   mask is .true. at gridpoints where band calculations are needed
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
!  <OUT NAME="ssalbband" TYPE="real">
!   the parameterization band values of the single-scattering albedo
!  </OUT>
!  <OUT NAME="asymmband" TYPE="real">
!   the parameterization band values of the asymmetry factor
!  </OUT>
! </SUBROUTINE>
!
subroutine esfsw_thickavg_1band (nband, nivl1, nivl2, nivls, nbands, &
                                 extivl, ssalbivl, asymmivl, solflxivl, &
                                 mask, extband, ssalbband, asymmband)

!---------------------------------------------------------------------
!    esfsw_thickavg_1band uses the thick-averaging technique to define the 
!    single-scattering properties of the specified parameterization band
!    spectral interval from the  specified spectral intervals of the 
!    particular scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer,                  intent(in)       :: nband
integer,                  intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl
real, dimension(:,:,:  ), intent(inout)      :: extband, ssalbband,   &
                                              asymmband
logical, dimension(:,:,:), intent(in)      :: mask

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nband       the sw parameterization band for which the optical
!                properties are being calculated
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    nivls       number of specified scattering spectral intervals      
!    nbands
!    extivl      specified spectral values of the extinction coefficient
!    asymmivl    the specified spectral values of the asymmetry     
!                factor                                           
!    solflxivl   the solar flux in each specified scattering spectral
!                interval                                         
!    mask        logical indicating the points at which the band values 
!                should be calculated       
!
!  intent(inout) variables:
!
!    ssalbivl    the specified spectral values of the single-       
!                scattering albedo                                   
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    ssalbband   the parameterization band values of the single-   
!                scattering albedo                                  
!    asymmband   the parameterization band values of the asymmetry   
!                factor                                               
!    
!--------------------------------------------------------------------
!  module variable:
!
!    solflxbandref_local  the solar flux in each parameterization band  
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_bands_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

      if (nbands /= size(solflxbandref_local,1)) then
        call error_mesg ('esfsw_bands_mod',  &
         'nbands argument not equal to number of SW bands registered', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
! call thickavg_1band with correct solar band

      call thickavg (nband, nivl1, nivl2, nivls, nbands, extivl, &
                     ssalbivl, asymmivl, solflxivl, &
                     solflxbandref_local(nband),    &
                     mask, extband, ssalbband, asymmband)
!--------------------------------------------------------------------

end subroutine esfsw_thickavg_1band

!####################################################################
! <SUBROUTINE NAME="esfsw_thinavg">
!  <OVERVIEW>
!   Subroutine to use thin-averaging technique to define band interval
!   single scattering properties.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thin-averaging technique to define the single-scattering    
! properties of the parameterization band spectral intervals from the  
! specified spectral intervals of the particular scatterer.            
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine esfsw_thinavg (nivl1    , nivl2     , nivls   ,   &
!                        nbands, &
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        extband  , ssalbband , asymmband)
!  </TEMPLATE>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="nivls" TYPE="integer">
!   number of specified scattering spectral intervals
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <IN NAME="nbands" TYPE="integer">
!   number of spectral bands
!  </IN>
!  <INOUT NAME="ssalbivl" TYPE="real">
!   the specified spectral values of the single-scattering albedo
!  </INOUT>
!  <IN NAME="asymmivl" TYPE="real">
!   the specified spectral values of the asymmetry factor
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
!  <OUT NAME="ssalbband" TYPE="real">
!   the parameterization band values of the single-scattering albedo
!  </OUT>
!  <OUT NAME="asymmband" TYPE="real">
!   the parameterization band values of the asymmetry factor
!  </OUT>
! </SUBROUTINE>
!
subroutine esfsw_thinavg (nivl1, nivl2, nivls, nbands, extivl, &
                          ssalbivl, asymmivl,  solflxivl, extband,   &
                          ssalbband , asymmband)

!---------------------------------------------------------------------
!    esfsw_thinavg uses the thin-averaging technique to define the 
!    single-scattering properties of the parameterization band spectral
!    intervals from the  specified spectral intervals of the particular
!    scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl
real, dimension(:,:,:,:), intent(out)      :: extband, ssalbband,   &
                                              asymmband

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    nivls       number of specified scattering spectral intervals      
!    nbands
!    extivl      specified spectral values of the extinction coefficient
!    asymmivl    the specified spectral values of the asymmetry     
!                factor                                           
!    solflxivl   the solar flux in each specified scattering spectral
!                interval                                         
!
!  intent(inout) variables:
!
!    ssalbivl    the specified spectral values of the single-       
!                scattering albedo                                   
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    ssalbband   the parameterization band values of the single-   
!                scattering albedo                                  
!    asymmband   the parameterization band values of the asymmetry   
!                factor                                               
!    
!--------------------------------------------------------------------
!  module variables:
!
!    solflxbandref_local  the solar flux in each parameterization band  
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_bands_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
! call thinavg with solar flux

      call thinavg (nivl1, nivl2, nivls, nbands, extivl, ssalbivl, &
                    asymmivl, solflxivl, solflxbandref_local,      &
                    extband, ssalbband, asymmband)
!---------------------------------------------------------------------

end subroutine esfsw_thinavg

!#####################################################################

  subroutine esfsw_bands_end


!--------------------------------------------------------------------
!    this is the destructor for esfsw_bands_mod
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('esfsw_bands_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    deallocate module data
!---------------------------------------------------------------------
      deallocate(endwvnbands_local)
      deallocate(solflxbandref_local)
      deallocate(solarfluxtoa_local)

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!-------------------------------------------------------------------

  end subroutine esfsw_bands_end

!#####################################################################

end module esfsw_bands_mod

