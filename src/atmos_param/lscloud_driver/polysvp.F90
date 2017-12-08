                         MODULE polysvp_mod

!------------------------------------------------------------------------
!
!         Contains interface routine compute_qs_a which calculates 
!              saturation vapor pressure and saturation mixing ratio 
!              either by:
!
!              a) calling FMS routine compute_qs;
!              b) using Goff-Gratch relations  contained within this
!                 module to define es with respect to water and ice, and 
!                 then using the standard (non-exact) expression for qs with
!                 respect to water and ice.
!
!-------------------------------------------------------------------------

use fms_mod,               only : FATAL, error_mesg, write_version_number
use sat_vapor_pres_mod,    only : compute_qs, sat_vapor_pres_init
use constants_mod,         only : RDGAS, RVGAS, HLV, HLS, TFREEZE, CP_AIR
use lscloud_constants_mod, only : lscloud_constants_init, mg_pr
use lscloud_types_mod,     only : lscloud_types_init, atmos_state_type,   &
                                  lscloud_nml_type

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------
!---interfaces-----------------------------------------------------------

interface polysvp_l
   module procedure polysvp_l1d, polysvp_l3d
end interface polysvp_l

interface polysvp_i
   module procedure polysvp_i1d, polysvp_i3d
end interface polysvp_i

PUBLIC  polysvp_init, compute_qs_a,  compute_qs_x1, &
        polysvp_l, polysvp_i, polysvp_end

!-------------------------------------------------------------------------
!----version number-------------------------------------------------------

Character(len=128) :: Version = '$Id$'
Character(len=128) :: Tagname = '$Name$'


!---------------------------------------------------------------------
real, parameter :: d622 = RDGAS / RVGAS
real, parameter :: d378 = 1. - d622

integer   :: super_ice_opt
logical   :: do_pdf_clouds
logical   :: pdf_org
logical   :: module_is_initialized = .false.



                              CONTAINS



!#########################################################################

SUBROUTINE polysvp_init (Nml_lsc)

type(lscloud_nml_type), intent(in) :: Nml_lsc

!-------------------------------------------------------------------------

      IF (module_is_initialized) return

      super_ice_opt = Nml_lsc%super_ice_opt
      do_pdf_clouds = Nml_lsc%do_pdf_clouds
      pdf_org = Nml_lsc%pdf_org

!------------------------------------------------------------------------
!    write version number to output file.
!------------------------------------------------------------------------
      call write_version_number (version, tagname)

!------------------------------------------------------------------------
!    make sure needed modules have been initialized.
!------------------------------------------------------------------------
      CALL sat_vapor_pres_init
      call lscloud_types_init
      call lscloud_constants_init

!------------------------------------------------------------------------
!    mark this module as initialized.
!------------------------------------------------------------------------
      module_is_initialized = .TRUE.


END SUBROUTINE polysvp_init



!#########################################################################

SUBROUTINE compute_qs_a     &
                (idim, jdim, kdim, tin, qin, pfull, ql_in, qi_in,   &
                                               conv_hum_ratio, Atmos_state)

!-------------------------------------------------------------------------
INTEGER,                INTENT(IN )   :: idim, jdim, kdim
real, dimension(:,:,:), intent(in)    :: tin, qin, pfull, ql_in, qi_in, &
                                         conv_hum_ratio
type(atmos_state_type), intent(inout) :: Atmos_state

!-----------------------------------------------------------------------
!----local variables

      INTEGER :: i, j, k

!-------------------------------------------------------------------------
!    compute qs and associated parameters, using different algorithms as
!    appropriate.    
!-------------------------------------------------------------------------
      if (do_pdf_clouds .and. pdf_org) then
        if (super_ice_opt .LT. 1) then

!-------------------------------------------------------------------------
!    if using pdf clouds and not allowing supersaturation with respect 
!    to ice:
!-------------------------------------------------------------------------
          call compute_qs (tin - ((HLV*ql_in + HLS*qi_in)/CP_AIR), pfull, &
                           Atmos_state%qs, dqsdT=Atmos_state%dqsdT,   &
                                                   esat=Atmos_state%esat0)
          Atmos_state%gamma = Atmos_state%dqsdT*   &
            (min(1., max(0., 0.05*(tin - ((HLV*ql_in +    &
                                                HLS*qi_in)/CP_AIR) -   &
                           TFREEZE + 20.)))*HLV +     &
             min(1., max(0., 0.05*(TFREEZE - tin + ((HLV*ql_in +   &
                                         HLS*qi_in)/CP_AIR))))*   &
                                                              HLS)/CP_AIR
        else

!-------------------------------------------------------------------------
!    if using pdf clouds and allowing supersaturation with respect to ice:
!-------------------------------------------------------------------------
          call error_mesg ( 'polysvp_mod/compute_qs_a', &
           'super_ice_opt > 1 when pdf_opt currently not supported', FATAL)
        end if
      else

!-------------------------------------------------------------------------
!    if not using pdf clouds and not allowing supersaturation with respect
!    to ice:
!-------------------------------------------------------------------------
        if (super_ice_opt .LT. 1) then

!------------------------------------------------------------------------
!    Calculate saturation specific humidity and its temperature 
!    derivative, and thermal conductivity plus vapor diffusivity factor.
!
!    FOR ORIGINAL SCHEME These are calculated according to the formulas:
!
!    (1)  qs   = d622*esat/ [pfull  -  (1.-d622)*esat]
!
!    (2) dqsdT = d622*pfull*(desat/dT)/[pfull-(1.-d622)*esat]**2.
!
!    (3) gamma = (L/cp) * dqsdT
!       
!       where d622 = rdgas/rvgas; esat = saturation vapor pressure;
!       and desat/dT is the temperature derivative of esat.
!       Note that in the calculation of gamma, 
!
!            {             hlv          for T > tfreeze             }
!       L =  { 0.05*(T-tfreeze+20.)*hlv + 0.05*(tfreeze-T)*hls      }
!            {                          for tfreeze-20.< T < tfreeze}
!            {             hls          for T < tfreeze-20.         }
!
!       This linear form is chosen because at tfreeze-20. es = esi, and
!       at tfreeze, es = esl, with linear interpolation in between.
!
!
!    (5) chi = 2.21 E-05 (m*m)/s  * (1.E+05)/pfull
!
!        where p is the pressure in Pascals.
!    
!
!       Note that qs, dqsdT, and gamma do not have their proper values
!       until all of the following code has been executed.  That
!       is qs and dqsdT are used to store intermediary results
!       in forming the full solution.
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!    calculate water saturated vapor pressure from table and store 
!    temporarily in the variable esat0.
!------------------------------------------------------------------------
          call compute_qs   &
                    (tin, pfull, Atmos_state%qs, dqsdT=Atmos_state%dqsdT, &
                                                    esat=Atmos_state%esat0)
          Atmos_state%gamma = Atmos_state%dqsdT*(min(1.,   &
            max(0., 0.05*(tin - TFREEZE + 20.)))*HLV +     &
             min(1., max(0., 0.05*(TFREEZE - tin)))*HLS)/CP_AIR
        else

!-------------------------------------------------------------------------
!    if not using pdf clouds and allowing supersaturation with respect 
!    to ice. in this case we need qs wrt ice and qs wrt liquid, which are
!    obtained in compute_qs_x1, albeit from a different formulation 
!    (Goff - Gratch) than is used elsewhere in FMS.  Consistency issues ???
!-------------------------------------------------------------------------
          call compute_qs_x1     &
               (idim, jdim, kdim, tin, pfull, Atmos_state%qs, &
                      Atmos_state%qsl, Atmos_state%qsi,   &
                                      Atmos_state%dqsdT, Atmos_state%gamma)
        end if
      end if

!--------------------------------------------------------------------------
!    Calculate relative humidity outside of convective cloud portion of 
!    grid box (U_ca) under the assumption that the temp is uniform across 
!    the gridbox. Define the maximum portion of the remaining grid-box
!    area which can be cloudy while maintaining the grid-box-mean relative
!    humidity (U01). Here conv_hum_ratio is the ratio of the gridbox RH to 
!    that in the environment of the convective cloud. 
!--------------------------------------------------------------------------
      IF (super_ice_opt .LT. 1) THEN

!--------------------------------------------------------------------------
!    when supersaturation not permitted, upper limit of 100% placed on 
!    external relative humidity.
!--------------------------------------------------------------------------
        where (conv_hum_ratio .gt. 0.)
          Atmos_state%U_ca =     &
               min(max(0., (qin/(conv_hum_ratio*Atmos_state%qs))), 1.)
        elsewhere
          Atmos_state%U_ca = 0.
        end where
        Atmos_state%U01 = Atmos_state%U_ca
      ELSE

!------------------------------------------------------------------------
!    when supersaturation is permitted, no upper limit placed on external
!    relative humidity. Max cloud fraction  depends on whether clouds are
!    either liquid or ice because of thedifference between qsl and qsi.
!--------------------------------------------------------------------------
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (conv_hum_ratio(i,j,k) .gt. 0. ) then
                Atmos_state%U_ca(i,j,k) = max(0.,  &
                     (qin(i,j,k)/(conv_hum_ratio(i,j,k)*   &
                                                  Atmos_state%qs(i,j,k))))
                if (tin(i,j,k) .LT. TFREEZE) then
                  Atmos_state%U01(i,j,k) = max(0.,  &
                     (qin(i,j,k)/(conv_hum_ratio(i,j,k)* &
                                                 Atmos_state%qsi(i,j,k))))
                else
                  Atmos_state%U01(i,j,k) = max(0.,  &
                     (qin(i,j,k)/(conv_hum_ratio(i,j,k)* &
                                                  Atmos_state%qsl(i,j,k))))
                endif
              else
                Atmos_state%U_ca(i,j,k) = 0.
              endif     
            end do
          end do
        end do
      END IF

!-----------------------------------------------------------------------


END SUBROUTINE compute_qs_a


!#########################################################################

SUBROUTINE compute_qs_x1 (idim, jdim, kdim, ttmp, pfull, qs, qs_l, qs_i,&
                          dqsdT, gamma )

!--------------------------------------------------------------------------
!    subroutine compute_qs_x1 computes es and qs over ice and liquid, using
!    the Goff-Gratch and Clausius-Clapeyron equations. Also calculated is
!    dqsdT and the term L/Cp * dqsdt (gamma).
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
integer,                         intent(in)     :: idim, jdim, kdim
REAL, dimension(idim,jdim,kdim), INTENT(IN )    :: ttmp, pfull
REAL, dimension(idim,jdim,kdim), INTENT(INOUT ) :: qs, dqsdT, gamma,  &
                                                   qs_l, qs_i

!--------------------------------------------------------------------------
!---local variables--------------------------------------------------------
      REAL, dimension (idim,jdim,kdim) :: eslt, esit
      integer                          :: i,j,k

!-------------------------------------------------------------------------
!    compute es with respect to liquid and ice.
!-------------------------------------------------------------------------
      eslt = polysvp_l3d (ttmp, idim, jdim, kdim)
      esit = polysvp_i3d (ttmp, idim, jdim, kdim)

!-------------------------------------------------------------------------
!    make sure ice saturation doesn't exceed water saturation at temps 
!    near freezing.
!-------------------------------------------------------------------------
      where (esit .gt. eslt) esit = eslt

!-------------------------------------------------------------------------
!    compute saturation specific humidity over liquid. calculate denominator
!    in qsat formula. limit denominator to esat, and thus qs to d622. this 
!    is done to avoid blow up in the upper stratosphere where pfull ~ esat.
!-------------------------------------------------------------------------
      qs_l = pfull - d378*eslt
      qs_l = max(qs_l, eslt) 
      qs_l = d622*eslt/qs_l

!-------------------------------------------------------------------------
!    compute saturation specific humidity over ice. calculate denominator
!    in qsat formula. limit denominator to esat, and thus qs to d622. this 
!    is done to avoid blow up in the upper stratosphere where pfull ~ esat.
!    compute saturation specific humidity over liquid.
!-------------------------------------------------------------------------
      qs_i = pfull - d378*esit
      qs_i = max(qs_i, esit) 
      qs_i = d622*esit/qs_i

!--------------------------------------------------------------------------
!    define the appropriate qs, dqsdT and gamma to use, dependent on ambient
!    temperature.
!--------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            IF (ttmp(i,j,k) .GE. TFREEZE - 23.16) THEN
              qs(i,j,k) = qs_l(i,j,k)
            ELSE
              qs(i,j,k) = qs_i(i,j,k)
            END IF
            IF (ttmp(i,j,k) .GE. TFREEZE - 40.) THEN
              dqsdT(i,j,k) = HLV*qs_l(i,j,k)/(RVGAS*ttmp(i,j,k)**2)
              gamma(i,j,k) = dqsdT(i,j,k)*HLV/CP_AIR
            ELSE
              dqsdT(i,j,k) = HLS*qs_i(i,j,k)/(RVGAS*ttmp(i,j,k)**2)
              gamma(i,j,k) = dqsdT(i,j,k)*HLS/CP_AIR
            END IF
          end do
        end do
     end do

!------------------------------------------------------------------------

END SUBROUTINE compute_qs_x1



!#########################################################################

function polysvp_l1d (T)

!------------------------------------------------------------------------
!    Compute saturation vapor pressure with respect to liquid  by using 
!    function from Goff and Gatch (1946). polysvp returned in units of pa.
!    T is input in units of K.
!------------------------------------------------------------------------

REAL(kind=mg_pr) :: T,polysvp_l1d

!-------------------------------------------------------------------------
!    Goff Gatch equation, uncertain below -70 C
!-------------------------------------------------------------------------
      polysvp_l1d = 10._mg_pr**(-7.90298_mg_pr*(373.16_mg_pr/t-1._mg_pr) + &
                    5.02808_mg_pr*log10(373.16_mg_pr/t) - &
                    1.3816e-7_mg_pr*(10._mg_pr**(11.344_mg_pr*(1._mg_pr - &
                                          t/373.16_mg_pr)) - 1._mg_pr) + &
                    8.1328e-3_mg_pr*(10._mg_pr**(-3.49149_mg_pr*   &
                            (373.16_mg_pr/t - 1._mg_pr)) - 1._mg_pr) + &
                                          log10(1013.246_mg_pr))*100._mg_pr 

!-------------------------------------------------------------------------


end function polysvp_l1d



!##########################################################################

function polysvp_i1d (T)

!------------------------------------------------------------------------
!    Compute saturation vapor pressure with respect to ice by using 
!    function from Goff and Gatch (1946). Polysvp returned in units of pa.
!    T is input in units of K.
!------------------------------------------------------------------------

REAL(kind=mg_pr) ::  T,polysvp_i1d

!-----------------------------------------------------------------------
! Goff Gatch equation  for ice (good down to -100 C)
!-----------------------------------------------------------------------
      polysvp_i1d = 10._mg_pr**(-9.09718_mg_pr*(273.16_mg_pr/t -   &
                                            1._mg_pr) - 3.56654_mg_pr* &
                     log10(273.16_mg_pr/t) + 0.876793_mg_pr*  &
                                         (1._mg_pr - t/273.16_mg_pr) + &
                                            log10(6.1071_mg_pr))*100._mg_pr

!-------------------------------------------------------------------------

end function polysvp_i1d


!########################################################################

function polysvp_l3d (T, idim, jdim, kdim)

!------------------------------------------------------------------------
!    Compute saturation vapor pressure with respect to liquid  by using 
!    function from Goff and Gatch (1946). polysvp returned in units of pa.
!    T is input in units of K.
!------------------------------------------------------------------------

REAL(kind=mg_pr), dimension(idim,jdim,kdim) :: T,polysvp_l3d
integer                                     :: idim,jdim,kdim

!-------------------------------------------------------------------------
!    Goff Gatch equation, uncertain below -70 C
!-------------------------------------------------------------------------
      polysvp_l3d = 10._mg_pr**(-7.90298_mg_pr*(373.16_mg_pr/t-1._mg_pr) + &
                    5.02808_mg_pr*log10(373.16_mg_pr/t) - &
                    1.3816e-7_mg_pr*(10._mg_pr**(11.344_mg_pr*(1._mg_pr - &
                                          t/373.16_mg_pr)) - 1._mg_pr) + &
                    8.1328e-3_mg_pr*(10._mg_pr**(-3.49149_mg_pr*   &
                            (373.16_mg_pr/t - 1._mg_pr)) - 1._mg_pr) + &
                                          log10(1013.246_mg_pr))*100._mg_pr 

!-------------------------------------------------------------------------


end function polysvp_l3d


!#########################################################################

function polysvp_i3d (T, idim, jdim, kdim)

!------------------------------------------------------------------------
!    Compute saturation vapor pressure with respect to ice by using 
!    function from Goff and Gatch (1946). Polysvp returned in units of pa.
!    T is input in units of K.
!------------------------------------------------------------------------

REAL(kind=mg_pr), dimension(idim,jdim,kdim) :: T, polysvp_i3d
integer                                     :: idim, jdim, kdim

!-----------------------------------------------------------------------
! Goff Gatch equation  for ice (good down to -100 C)
!-----------------------------------------------------------------------
      polysvp_i3d = 10._mg_pr**(-9.09718_mg_pr*(273.16_mg_pr/t -   &
                                            1._mg_pr) - 3.56654_mg_pr* &
                     log10(273.16_mg_pr/t) + 0.876793_mg_pr*  &
                                         (1._mg_pr - t/273.16_mg_pr) + &
                                            log10(6.1071_mg_pr))*100._mg_pr

!-------------------------------------------------------------------------
   
end function polysvp_i3d


!#########################################################################

SUBROUTINE  polysvp_end

      module_is_initialized = .FALSE.

END SUBROUTINE  polysvp_end


!#########################################################################






                     END MODULE polysvp_mod
