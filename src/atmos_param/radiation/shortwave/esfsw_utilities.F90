 
module esfsw_utilities_mod


use fms_mod,      only: write_version_number, &
                        error_mesg, FATAL
!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!-------  version number --------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------

public :: esfsw_utilities_init, &
          esfsw_utilities_end, &
          thickavg, thinavg

interface thickavg
   module procedure thickavg_3d
   module procedure thickavg_0d
   module procedure thickavg_1band
   module procedure thickavg_isccp
end interface

!---------------------------------------------------------------------
!------- private data ------

logical :: module_is_initialized=.false.   ! module is initialized ?

!---------------------------------------------------------------------

CONTAINS

!#####################################################################
!
!                     PUBLIC SUBROUTINES
!
!#####################################################################

subroutine esfsw_utilities_init

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    write version number to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)

!-------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------

end subroutine esfsw_utilities_init

!#####################################################################

subroutine esfsw_utilities_end

!--------------------------------------------------------------------
!    this is the destructor for esfsw_utilities_mod
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('esfsw_utilites_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!-------------------------------------------------------------------

end subroutine esfsw_utilities_end

!#####################################################################

! <SUBROUTINE NAME="thickavg_3d">
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
!   call subroutine thickavg_3d (nivl1    , nivl2     , nivls   ,   &
!                        nbands, $
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        solflxband, extband  , ssalbband , asymmband)
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
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
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
subroutine thickavg_3d (nivl1, nivl2, nivls, nbands, extivl, ssalbivl,&
                        asymmivl, solflxivl, solflxband, mask, extband,&
                        ssalbband, asymmband)
 
!---------------------------------------------------------------------
!    thickavg_3d uses the thick-averaging technique to define the 
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
real, dimension(:),       intent(in)       :: solflxband            
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
!    solflxband  the solar flux in each parameterization band  
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
!  local variables:
 
      real, dimension (size(ssalbivl,1),   &
                       size(ssalbivl,2), &
                       size(ssalbivl,3), nbands)  ::   refband        

      real, dimension (size(ssalbivl,1),   &
                       size(ssalbivl,2), &
                       size(ssalbivl,3))  ::   refthick, sp, sumk,   &
                                               sumomegak, sumomegakg, &
                                               sumrefthick

      integer  :: nband
      integer  :: i, j, k, ni
 
!--------------------------------------------------------------------
!  local variables:
!
!     refband
!     refthick
!     sp
!     sumk
!     sumomegak
!     sumomegakg
!     sumrefthck
!     nband
!     i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!------------------------------------------------------ --------------
!--------------------------------------------------------------------
      do nband = 1,nbands
        sumk(:,:,:) = 0.0
        sumomegak(:,:,:) = 0.0
        sumomegakg(:,:,:) = 0.0
        sumrefthick(:,:,:) = 0.0
        do ni = nivl1(nband),nivl2(nband)
!
          do k=1, size(ssalbivl,3)
            do j=1,size(ssalbivl,2)
              do i=1,size(ssalbivl,1)
                if (mask(i,j,k)) then
                  ssalbivl(i,j,k,ni) = MIN(ssalbivl(i,j,k,ni), 1.0)
                  sp(i,j,k) = sqrt( ( 1.0 - ssalbivl(i,j,k,ni) ) /    &
                                    ( 1.0 - ssalbivl(i,j,k,ni) *      &
                                      asymmivl(i,j,k,ni) ) )
                  refthick(i,j,k) = (1.0 - sp(i,j,k))/(1.0 + sp(i,j,k))
                  sumrefthick(i,j,k) = sumrefthick(i,j,k) +    &
                                       refthick(i,j,k)*  &
                                       solflxivl(nband,ni)
                  sumk(i,j,k) = sumk(i,j,k) + extivl(i,j,k,ni) *   &
                                solflxivl(nband,ni)
                  sumomegak(i,j,k) = sumomegak(i,j,k) +     &
                                     ssalbivl(i,j,k,ni)*   &
                                     extivl(i,j,k,ni) *   &
                                     solflxivl(nband,ni)
                  sumomegakg(i,j,k) = sumomegakg(i,j,k) +    &
                                      ssalbivl(i,j,k,ni)*&
                                      extivl(i,j,k,ni)*  &
                                      asymmivl(i,j,k,ni) * &
                                      solflxivl(nband,ni)
                endif
              end do
            end do
          end do
        end do

!---------------------------------------------------------------------
!    the 1.0E-100 factor to calculate asymmband is to prevent        
!    division by zero.                                             
!---------------------------------------------------------------------
        do k=1, size(ssalbivl,3)
          do j=1,size(ssalbivl,2)
            do i=1,size(ssalbivl,1)
              extband(i,j,k,nband) = sumk(i,j,k) / solflxband(nband)
              asymmband(i,j,k,nband) = sumomegakg(i,j,k) /         &
                                       ( sumomegak(i,j,k) + 1.0E-100)
              refband(i,j,k,nband) = sumrefthick(i,j,k)/  &
                                     solflxband(nband)
              ssalbband(i,j,k,nband) = 4.0 * refband(i,j,k,nband) / &
                                       ((1.0 +    &
                                       refband(i,j,k,nband)) ** 2 -&
                                       asymmband(i,j,k,nband) *     &
                                       (1.0 - refband(i,j,k,nband))**2 )
            end do
          end do
        end do
      end do

!---------------------------------------------------------------------

  
end subroutine thickavg_3d



!####################################################################
! <SUBROUTINE NAME="thickavg_0d">
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
!   call subroutine thickavg_0d (nivl1    , nivl2     , nivls   ,   &
!                        nbands,  &
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        solflxband, extband  , ssalbband , asymmband)
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
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
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
subroutine thickavg_0d (nivl1, nivl2, nivls, nbands, extivl, ssalbivl,&
                        asymmivl, solflxivl, solflxband, extband,  &
                        ssalbband , asymmband)
 
!---------------------------------------------------------------------
!    thickavg_0d uses the thick-averaging technique to define the 
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
real, dimension(:),       intent(in)       :: solflxband            
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
!    solflxband  the solar flux in each parameterization band  
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
!  local variables:

      real, dimension (nbands)        ::   refband
      real                            ::   refthick, sp, sumk,   &
                                           sumomegak, sumomegakg, &
                                           sumrefthick
      integer  :: nband 
      integer  :: ni
 
 
!--------------------------------------------------------------------
!  local variables:
!
!     refband
!     refthick
!     sp
!     sumk
!     sumomegak
!     sumomegakg
!     sumrefthck
!     nband
!     i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif
  
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      do nband = 1,nbands
        sumk        = 0.0
        sumomegak   = 0.0
        sumomegakg  = 0.0
        sumrefthick = 0.0
        do ni = nivl1(nband),nivl2(nband)
          if (extivl(ni) /= 0.0) then
            ssalbivl(ni) = MIN(ssalbivl(ni), 1.0)
            sp = sqrt( ( 1.0 - ssalbivl(ni) ) /    &
                       ( 1.0 - ssalbivl(ni) * asymmivl(ni) ) )
            refthick = (1.0 - sp)/(1.0 + sp)
            sumrefthick = sumrefthick + refthick * solflxivl(nband,ni)
            sumk = sumk + extivl(ni) * solflxivl(nband,ni)
            sumomegak = sumomegak +     &
                        ssalbivl(ni) * extivl(ni) * solflxivl(nband,ni)
            sumomegakg = sumomegakg +    &
                         ssalbivl(ni) * extivl(ni) *  &
                         asymmivl(ni) * solflxivl(nband,ni)
          endif
        end do

!--------------------------------------------------------------------- 
!
!--------------------------------------------------------------------- 
        extband(nband) = sumk / solflxband(nband)
        asymmband(nband) = sumomegakg / ( sumomegak + 1.0E-100)
        refband(nband) = sumrefthick/ solflxband(nband)
        ssalbband(nband) = 4.0 * refband(nband) / &
                           ( (1.0 + refband(nband))**2 - &
                          asymmband(nband) * (1.0 - refband(nband))**2 )
      end do

!---------------------------------------------------------------------
  

end subroutine thickavg_0d


!##################################################################

! <SUBROUTINE NAME="thickavg_isccp">
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
!   call subroutine thickavg (nband, nivl1, nivl2, extivl, solflxivl, &
!                             solflxband, mask, extband )
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
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
!  </IN>
!  <IN NAME="mask" TYPE="logical">
!   mask is .true. at gridpoints where extband needs to be calculated
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
! </SUBROUTINE>
!
subroutine thickavg_isccp (nband, nivl1, nivl2, extivl,          &
                           solflxivl, solflxband, mask, extband)
 
!---------------------------------------------------------------------
!    thickavg_isccp uses the thick-averaging technique to define the 
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
real,                     intent(in)       :: solflxband            
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
!    solflxband  the solar flux in each parameterization band  
!    mask        logical indicating the points at which the band values 
!                should be calculated       
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!  local variables:
 
      real     ::  sumk
      integer  ::  i, j, k, ni
 
!--------------------------------------------------------------------
!  local variables:
!
!     sumk
!     i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!
      do k=1, size(extivl,3)
        do j=1,size(extivl,2)
          do i=1,size(extivl,1)
            if (mask(i,j,k)) then
              sumk = 0.0
              do ni = nivl1,nivl2
                sumk = sumk + extivl(i,j,k,ni)*solflxivl(nband,ni)
              end do
              extband(i,j,k) = sumk/solflxband
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------

  
end subroutine thickavg_isccp

!##################################################################

! <SUBROUTINE NAME="thickavg_1band">
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
!   call subroutine thickavg (nband, nivl1  , nivl2   , nivls   , &
!                        nbands, &
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        solflxband,  mask, extband  , ssalbband ,  &
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
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
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
subroutine thickavg_1band (nband, nivl1, nivl2, nivls, nbands, extivl, &
                           ssalbivl, asymmivl, solflxivl, solflxband, &
                           mask, extband, ssalbband, asymmband)
 
!---------------------------------------------------------------------
!    thickavg_1band uses the thick-averaging technique to define the 
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
real,                     intent(in)       :: solflxband            
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
!    solflxband  the solar flux in each parameterization band  
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
 
!--------------------------------------------------------------------
!  local variables:
 
      real :: refband, sp, refthick
      real :: sumk, sumomegak, sumomegakg,  sumrefthick

      integer  :: i, j, k, ni
 
!--------------------------------------------------------------------
!  local variables:
!
!     refband
!     refthick
!     sp
!     sumk
!     sumomegak
!     sumomegakg
!     sumrefthck
!     nband
!     i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      do k=1, size(ssalbivl,3)
        do j=1,size(ssalbivl,2)
          do i=1,size(ssalbivl,1)
            if (mask(i,j,k)) then
              sumk        = 0.0
              sumomegak        = 0.0
              sumomegakg        = 0.0
              sumrefthick        = 0.0
              do ni = nivl1,nivl2
                ssalbivl(i,j,k,ni) = MIN(ssalbivl(i,j,k,ni), 1.0)
                sp = sqrt((1.0 - ssalbivl(i,j,k,ni) ) /    &
                          (1.0 - ssalbivl(i,j,k,ni)*asymmivl(i,j,k,ni)))
                refthick = (1.0 - sp)/(1.0 + sp)
                sumrefthick = sumrefthick + refthick*solflxivl(nband,ni)
                sumk = sumk + extivl(i,j,k,ni)*solflxivl(nband,ni)
                sumomegak = sumomegak + ssalbivl(i,j,k,ni)*   &
                                        extivl(i,j,k,ni)*   &
                                        solflxivl(nband,ni)
                sumomegakg = sumomegakg + ssalbivl(i,j,k,ni)*&
                                          extivl(i,j,k,ni)*  &
                                          asymmivl(i,j,k,ni)* &
                                          solflxivl(nband,ni)
              end do

!---------------------------------------------------------------------
!    the 1.0E-100 factor to calculate asymmband is to prevent        
!    division by zero.                                             
!---------------------------------------------------------------------
              extband(i,j,k) = sumk/solflxband
              asymmband(i,j,k) = sumomegakg/(sumomegak + 1.0E-100)
              refband  = sumrefthick/solflxband        
              ssalbband(i,j,k) = 4.0*refband/((1.0 + refband) ** 2 - &
                                 asymmband(i,j,k)*(1.0 - refband)**2 )
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------

  
end subroutine thickavg_1band

!####################################################################
! <SUBROUTINE NAME="thinavg">
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
!   call subroutine thinavg (nivl1    , nivl2     , nivls   ,   &
!                        nbands, &
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        solflxband, extband  , ssalbband , asymmband)
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
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
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
subroutine thinavg (nivl1, nivl2, nivls, nbands, extivl, ssalbivl, &
                    asymmivl,  solflxivl, solflxband, extband,   &
                    ssalbband , asymmband)
 
!---------------------------------------------------------------------
!    thinavg uses the thin-averaging technique to define the 
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
real, dimension(:),       intent(in)       :: solflxband            
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
!    solflxband  the solar flux in each parameterization band  
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
!  local variables:

      real, dimension (size(ssalbivl,1),   &
                       size(ssalbivl,2),  &
                       size(ssalbivl,3)) ::   sumk,  sumomegak,   &
                                              sumomegakg
 
      integer   ::   nband
      integer   ::   i, j, k, ni

!--------------------------------------------------------------------
!  local variables:
! 
!    sumk
!    sumomegak
!    sumomegakg
!    nband
!    i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('esfsw_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do nband = 1,nbands
        sumk(:,:,:) = 0.0
        sumomegak(:,:,:) = 0.0
        sumomegakg(:,:,:) = 0.0
        do ni = nivl1(nband),nivl2(nband)
          do k=1, size(ssalbivl,3)
            do j=1,size(ssalbivl,2)
              do i=1,size(ssalbivl,1)
                if ((ssalbivl(i,j,k,ni) +    &
                     asymmivl(i,j,k,ni)) /= 0.0) then
                  ssalbivl(i,j,k,ni) = MIN(ssalbivl(i,j,k,ni), 1.0)
                  sumk(i,j,k) = sumk(i,j,k) + extivl(i,j,k,ni) *   &
                                solflxivl(nband,ni)
                  sumomegak(i,j,k) = sumomegak(i,j,k) +    &
                                     ssalbivl(i,j,k,ni) *  &
                                     extivl(i,j,k,ni) *   &
                                     solflxivl(nband,ni)
                  sumomegakg(i,j,k) = sumomegakg(i,j,k) +    &
                                      ssalbivl(i,j,k,ni) * & 
                                      extivl(i,j,k,ni) *   &
                                      asymmivl(i,j,k,ni) *  &
                                      solflxivl(nband,ni)
                endif
              end do
            end do
          end do
        end do

!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
        do k=1, size(ssalbivl,3)
          do j=1,size(ssalbivl,2)
            do i=1,size(ssalbivl,1)
              extband(i,j,k,nband) = sumk(i,j,k) / solflxband(nband)
              asymmband(i,j,k,nband) = sumomegakg(i,j,k) /    &
                                       ( sumomegak(i,j,k) + 1.0E-100 )
              ssalbband(i,j,k,nband) = sumomegak(i,j,k) /   &
                                       ( sumk(i,j,k) + 1.0E-100 )
            end do
          end do
        end do
      end do

!-------------------------------------------------------------------
  

end subroutine thinavg 


!#####################################################################

end module esfsw_utilities_mod

