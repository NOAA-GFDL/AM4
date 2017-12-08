module cloud_chem

  use constants_mod, only : AVOGNO, PSTD_MKS

  !references
  ! jacobson (2005)
  implicit none

  public :: cloud_ph, cloud_so2_chem, cloud_nb_diag

  public :: CLOUD_CHEM_PH_LEGACY, CLOUD_CHEM_PH_BISECTION, CLOUD_CHEM_PH_CUBIC
  public :: CLOUD_CHEM_LEGACY, CLOUD_CHEM_F1P, CLOUD_CHEM_F1P_BUG, CLOUD_CHEM_F1P_BUG2

  private

  real*8, parameter      ::  ra = 8314./PSTD_MKS   ! universal constant   (atm)/(m-k)
  integer                ::  n_h2o2, n_hno3, n_hcooh, n_ch3cooh, n_so2, n_nh3 !tracer indices
  integer, parameter     ::  cloud_nb_diag = 11

  integer, parameter     :: CLOUD_CHEM_PH_LEGACY    = 1
  integer, parameter     :: CLOUD_CHEM_PH_BISECTION = 2
  integer, parameter     :: CLOUD_CHEM_PH_CUBIC     = 3

  integer, parameter     :: CLOUD_CHEM_LEGACY    = 1
  integer, parameter     :: CLOUD_CHEM_F1P       = 2
  integer, parameter     :: CLOUD_CHEM_F1P_BUG   = 3
  integer, parameter     :: CLOUD_CHEM_F1P_BUG2  = 4
  
  real*8, parameter      :: const0 = 1.e3/AVOGNO

contains

  subroutine cloud_pH(ihno3,iso2,inh3,ihcooh,ich3cooh,iso4,ico2,ialk,tk,patm,xlw,air_dens,cloud_solver,ch,diagnostics)

    !input

    real*8, intent(in)     :: ihno3,iso2,inh3,ihcooh,ich3cooh !vmr
    real*8, intent(in)     :: iso4,ico2,ialk                  !vmr
    real*8, intent(in)     :: xlw                             !l(water)/l(air)
    real*8, intent(in)     :: air_dens                        !molec/cm3
    real*8, intent(in)     :: tk                              !t[k]
    real*8, intent(in)     :: patm                            !atmospheric pressure (atm)
    integer, intent(in)    :: cloud_solver                

    !ouput
    real*8, intent(out)    :: ch                         !h concentration in mol/l
    real*8, intent(out),optional, dimension(cloud_nb_diag)  :: diagnostics


    if ( cloud_solver .eq. cloud_chem_ph_legacy ) then
       call cloud_pH_am3(ihno3,iso2,inh3,ihcooh,ich3cooh,iso4,ico2,tk,patm,xlw,air_dens,ch,diagnostics)
    elseif ( cloud_solver .eq. cloud_chem_ph_cubic ) then
       call cloud_pH_cubic(ihno3,iso2,inh3,ihcooh,ich3cooh,iso4,ico2,tk,patm,xlw,air_dens,ch,diagnostics)
    elseif  ( cloud_solver .eq. cloud_chem_ph_bisection ) then
       call cloud_pH_bisection(ihno3,iso2,inh3,ihcooh,ich3cooh,iso4,ico2,ialk,tk,patm,xlw,air_dens,ch,diagnostics)
    end if

  end subroutine cloud_pH

  subroutine cloud_pH_am3(ihno3,iso2,inh3,ihcooh,ich3cooh,iso4,ico2,tk,patm,xlw,air_dens,ch,diagnostics)

    real*8, parameter      :: ch0 = 1.e-5       
    integer, parameter   :: iter_max = 40

    !input

    real*8, intent(in)     :: ihno3,iso2,inh3,ihcooh,ich3cooh ! vmr
    real*8, intent(in)     :: iso4,ico2                       ! vmr
    real*8, intent(in)     :: xlw                             ! l(water)/l(air)
    real*8, intent(in)     :: air_dens                        ! molec/cm3
    real*8, intent(in)     :: tk                              ! t[k]
    real*8, intent(in)     :: patm                            ! atmospheric pressure (atm)

    real*8, intent(out),optional, dimension(cloud_nb_diag)  :: diagnostics

    !ouput

    real*8, intent(out)    :: ch                              ! h concentration in mol/l

    !local variables
    real*8                 :: delta

    integer              :: iter
    logical              :: converged
    real*8                 :: kw
    real*8                 :: he_hno3,he_hcooh,he_ch3cooh,he_nh3, he_so2
    real*8                 :: h_hno3,h_hcooh,h_ch3cooh,h_nh3, h_so2, h_co2
    real*8                 :: ka_hno3,ka_hcooh,ka_ch3cooh,kb_nh3, ka1_so2, ka2_so2, ka1_co2,ka2_co2
    real*8                 :: h_x_ka_hno3,h_x_ka1_so2,h_x_ka1_x_ka2_so2,h_x_kb_kw_nh3,h_x_ka_hcooh,h_x_ka_ch3cooh,h_x_ka1_co2,h_x_ka1_x_ka2_co2
    real*8                 :: e_hno3, e_so2, e_so4, e_hcooh, e_ch3cooh, e_nh3, e_h2o, e_co2
    real*8                 :: hno3g, so2g, hcoohg, ch3coohg, nh3g, co2g
    real*8                 :: hno3, so2, hcooh, ch3cooh, nh3, co2, so4
    real*8                 :: conv, chm1,com2

    ch = ch0                      ! initial h concentration

    call so2aq(tk,h_so2,ka1_so2,ka2_so2)
    call hno3aq(tk,h_hno3,ka_hno3)
    call nh3aq(tk,h_nh3,kb_nh3)
    call hcoohaq(tk,h_hcooh,ka_hcooh)
    call ch3coohaq(tk,h_ch3cooh,ka_ch3cooh)
    call co2aq(tk,h_co2,ka1_co2,ka2_co2)
    call kw_aq(tk,kw)
    
    conv = ra*tk*xlw
    
    h_x_ka_hno3               = h_hno3*ka_hno3
    h_x_ka1_so2               = h_so2*ka1_so2
    h_x_ka1_x_ka2_so2         = h_so2*ka1_so2 * ka2_so2
    h_x_kb_kw_nh3             = h_nh3 * kb_nh3/kw
    h_x_ka_hcooh              = h_hcooh*ka_hcooh
    h_x_ka_ch3cooh            = h_ch3cooh * ka_ch3cooh
    h_x_ka1_co2               = h_co2*ka1_co2
    h_x_ka1_x_ka2_co2         = h_co2*ka1_co2*ka2_co2


    hno3    = max(ihno3,0.)
    so2     = max(iso2,0.)
    nh3     = max(inh3,0.)
    hcooh   = max(ihcooh,0.)
    ch3cooh = max(ich3cooh,0.)
    co2     = max(ico2,0.)
    so4     = max(iso4,0.)
    
    !independent of ph
    e_so4 = so4*air_dens               &                       ! /cm3(a)
         *const0/xlw
       
    iter  = 0
    
    do while ( iter .lt. iter_max )

!       write(*,*) 'ch',ch

       delta = ch
       
       chm1 = 1./ch
       
       !calculate effective henry's constant          
       he_hno3     = h_hno3       +h_x_ka_hno3    *chm1
       he_so2      = h_so2        +h_x_ka1_so2    *chm1      +h_x_ka1_x_ka2_so2*chm1*chm1
       he_nh3      = h_nh3        +h_x_kb_kw_nh3  *ch
       he_hcooh    = h_hcooh      +h_x_ka_hcooh   *chm1
       he_ch3cooh  = h_ch3cooh    +h_x_ka_ch3cooh *chm1
       
!       write(*,*) 'he_nh3',he_nh3
!       write(*,*) 'he_hno3',he_hno3

       !gas concentration in equilibrium with cloud water
       hno3g       = hno3    / (1. + he_hno3     * conv)
       so2g        = so2     / (1. + he_so2      * conv)
       nh3g        = nh3     / (1. + he_nh3      * conv)
       hcoohg      = hcooh   / (1. + he_hcooh    * conv)
       ch3coohg    = ch3cooh / (1. + he_ch3cooh  * conv)
       co2g        = co2
       
       !concentrations in aqueous phase
       e_hno3     = h_x_ka_hno3    * hno3g    * patm
       e_so2      = h_x_ka1_so2    * so2g     * patm
       e_nh3      = h_x_kb_kw_nh3  * nh3g     * patm
       e_h2o      = kw
       e_co2      = h_x_ka1_co2    * co2g     * patm
       e_hcooh    = h_x_ka_hcooh   * hcoohg   * patm              
       e_ch3cooh  = h_x_ka_ch3cooh * ch3coohg * patm

       com2 = (e_hno3 + e_hcooh + e_ch3cooh + e_so2 + e_h2o + e_co2)/(1.+e_nh3)

       com2 = max(com2,1.e-20)
       
       ch    = sqrt(com2)

       ch    = min(1.e-2,max(1.e-7,ch+2*e_so4))
       
       if ( iter > 1 ) then
          delta = abs( ( ch - delta ) /delta )
          converged = delta < 1e-5
          if ( converged ) then
             exit
          end if
       end if
       
       iter = iter + 1
       
    end do
    
    if( .not. converged ) then
       write(*,*) 'pH not converged  (% change)=', &
            100.*delta
    end if

    if ( present(diagnostics) ) then
       chm1 = 1./ch
       
       !calculate effective henry's constant          
       he_hno3     = h_hno3       +h_x_ka_hno3    *chm1
       he_so2      = h_so2        +h_x_ka1_so2    *chm1      +h_x_ka1_x_ka2_so2*chm1*chm1
       he_nh3      = h_nh3        +h_x_kb_kw_nh3  *ch
       he_hcooh    = h_hcooh      +h_x_ka_hcooh   *chm1
       he_ch3cooh  = h_ch3cooh    +h_x_ka_ch3cooh *chm1
       
!       write(*,*) 'he_nh3',he_nh3
!       write(*,*) 'he_hno3',he_hno3

       !gas concentration in equilibrium with cloud water
       hno3g       = hno3    / (1. + he_hno3     * conv)
       so2g        = so2     / (1. + he_so2      * conv)
       nh3g        = nh3     / (1. + he_nh3      * conv)
       hcoohg      = hcooh   / (1. + he_hcooh    * conv)
       ch3coohg    = ch3cooh / (1. + he_ch3cooh  * conv)
       co2g        = co2

       diagnostics(1) = h_x_ka_hno3    * hno3g    * patm / ch
       diagnostics(2) = h_x_ka_hcooh   * hcoohg   * patm / ch
       diagnostics(3) = h_x_ka_ch3cooh * ch3coohg  * patm / ch
       diagnostics(4) = h_x_ka1_so2    * so2g      *patm / ch
       diagnostics(5) = h_x_ka1_x_ka2_so2  * so2g  * patm / (ch*ch)
       diagnostics(6) = h_x_ka1_co2  * co2g        * patm / ch
       diagnostics(7) = h_x_ka1_x_ka2_co2  * co2g  * patm / (ch*ch)
       diagnostics(8) = e_so4
       diagnostics(9) = h_x_kb_kw_nh3  * nh3g     * patm * ch
       diagnostics(10) = kw / ch
       diagnostics(11) = sum(diagnostics(1:4)) + 2*diagnostics(5) + diagnostics(6) + 2*diagnostics(7)+2*diagnostics(8) - diagnostics(9) + diagnostics(10) - cH

    end if

        
  end subroutine cloud_pH_am3


  subroutine cloud_pH_cubic(ihno3,iso2,inh3,ihcooh,ich3cooh,iso4,ico2,tk,patm,xlw,air_dens,ch,diagnostics)

    real*8, parameter      :: ch0 = 1.e-5       
    integer, parameter     :: iter_max = 100

    !input

    real*8, intent(in)     :: ihno3,iso2,inh3,ihcooh,ich3cooh ! vmr
    real*8, intent(in)     :: iso4,ico2                       ! vmr
    real*8, intent(in)     :: xlw                             ! l(water)/l(air)
    real*8, intent(in)     :: air_dens                        ! molec/cm3
    real*8, intent(in)     :: tk                              ! t[k]
    real*8, intent(in)     :: patm                            ! atmospheric pressure (atm)

    real*8, intent(out),optional, dimension(cloud_nb_diag)  :: diagnostics

    !ouput

    real*8, intent(out)    :: ch                              ! h concentration in mol/l

    !local variables
    real*8                 :: delta

    integer              :: iter
    logical              :: converged
    real*8                 :: kw
    real*8                 :: he_hno3,he_hcooh,he_ch3cooh,he_nh3, he_so2, he_co2
    real*8                 :: h_hno3,h_hcooh,h_ch3cooh,h_nh3, h_so2, h_co2
    real*8                 :: ka_hno3,ka_hcooh,ka_ch3cooh,kb_nh3, ka1_so2, ka2_so2, ka1_co2, ka2_co2
    real*8                 :: h_x_ka_hno3,h_x_ka1_so2,h_x_ka1_x_ka2_so2,h_x_kb_kw_nh3,h_x_ka_hcooh,h_x_ka_ch3cooh,h_x_ka1_co2,h_x_ka1_x_ka2_co2
    real*8                 :: e_hno3, e_so2_1,e_so2_2, e_so4, e_hcooh, e_ch3cooh, e_nh3, e_h2o, e_co2_1,e_co2_2
    real*8                 :: hno3g, so2g, hcoohg, ch3coohg, nh3g, co2g
    real*8                 :: hno3, so2, hcooh, ch3cooh, nh3, co2, so4
    real*8                 :: conv, chm1,com2

    !cubic solver
    real*8                 :: a0,a1,a2
    real*8                 :: a,b,d,e,f,g,h,p,q,r,crutes(3)
    integer              :: nr

    ch = ch0                      ! initial h concentration

    call so2aq(tk,h_so2,ka1_so2,ka2_so2)
    call hno3aq(tk,h_hno3,ka_hno3)
    call nh3aq(tk,h_nh3,kb_nh3)
    call hcoohaq(tk,h_hcooh,ka_hcooh)
    call ch3coohaq(tk,h_ch3cooh,ka_ch3cooh)
    call co2aq(tk,h_co2,ka1_co2,ka2_co2)
    call kw_aq(tk,kw)
    
    conv = ra*tk*xlw
    
    h_x_ka_hno3               = h_hno3*ka_hno3
    h_x_ka1_so2               = h_so2*ka1_so2
    h_x_ka1_x_ka2_so2         = h_so2*ka1_so2 * ka2_so2
    h_x_kb_kw_nh3             = h_nh3 * kb_nh3/kw
    h_x_ka_hcooh              = h_hcooh*ka_hcooh
    h_x_ka_ch3cooh            = h_ch3cooh * ka_ch3cooh
    h_x_ka1_co2               = h_co2*ka1_co2
    h_x_ka1_x_ka2_co2         = h_co2*ka1_co2*ka2_co2

    hno3    = max(ihno3,0.)
    so2     = max(iso2,0.)
    nh3     = max(inh3,0.)
    hcooh   = max(ihcooh,0.)
    ch3cooh = max(ich3cooh,0.)
    co2     = max(ico2,0.)
    so4     = max(iso4,0.)
    
    !independent of ph
    e_so4 = so4*air_dens               &                       ! /cm3(a)
         *const0/xlw

    d     = 2*e_so4
       
    iter  = 1
    
    do while ( iter .lt. iter_max )

!       write(*,*) 'ch',ch

       delta = ch
       
       chm1 = 1./ch
       
       !calculate effective henry's constant          
       he_hno3     = h_hno3       +h_x_ka_hno3    *chm1
       he_so2      = h_so2        +h_x_ka1_so2    *chm1      +h_x_ka1_x_ka2_so2*chm1*chm1
       he_nh3      = h_nh3        +h_x_kb_kw_nh3  *ch
       he_hcooh    = h_hcooh      +h_x_ka_hcooh   *chm1
       he_ch3cooh  = h_ch3cooh    +h_x_ka_ch3cooh *chm1
       he_co2      = h_co2        +h_x_ka1_co2    *chm1      +h_x_ka1_x_ka2_co2*chm1*chm1
       
       !gas concentration in equilibrium with cloud water
       hno3g       = hno3    / (1. + he_hno3     * conv)
       so2g        = so2     / (1. + he_so2      * conv)
       nh3g        = nh3     / (1. + he_nh3      * conv)
       hcoohg      = hcooh   / (1. + he_hcooh    * conv)
       ch3coohg    = ch3cooh / (1. + he_ch3cooh  * conv)
       co2g        = co2     / (1. + he_co2      * conv)
       
       !concentrations in aqueous phase
       e_hno3     = h_x_ka_hno3    * hno3g    * patm
       e_so2_1    = h_x_ka1_so2    * so2g     * patm
       e_so2_2    = 2. * e_so2_1   * ka2_so2
       e_nh3      = h_x_kb_kw_nh3  * nh3g     * patm
       e_h2o      = kw
       e_co2_1    = h_x_ka1_co2     * co2g     * patm
       e_co2_2    = 2. * e_co2_1   * ka2_co2
       e_hcooh    = h_x_ka_hcooh   * hcoohg   * patm              
       e_ch3cooh  = h_x_ka_ch3cooh * ch3coohg * patm

       e  = e_h2o + e_co2_1 + e_so2_1 + e_hno3 + e_hcooh + e_ch3cooh
       f  = e_co2_2 + e_so2_2
       g  = e_nh3
       
       h  = 1. + g
       p  = - d / h
       q  = - e / h
       r  =   f / h
       a  =  1./3.  * (3.*q - p*p)
       b  =  1./27. * (2.*p*p*p-9.*p*q+27.*r)

       a2  = 0.
       a1  = a
       a0  = b

       call cubic(a2,a1,a0,nr,crutes)

       ch  = crutes(1) - p/3.
       
       ch    = min(1e-1,max(1.e-13,ch))
       
       if ( iter > 1 ) then
          delta = abs( ( ch - delta ) /delta )
          converged = delta < 1e-5
          if ( converged ) then
             exit
          end if
       end if
       
       iter = iter + 1
       
    end do
    
    if( .not. converged ) then
       write(*,*) 'setsox: ph failed to o converge % change=', &
            100.*delta, 'final ch',ch
    end if

    if ( present(diagnostics) ) then

       chm1 = 1./ch
       
       !calculate effective henry's constant          
       he_hno3     = h_hno3       +h_x_ka_hno3    *chm1
       he_so2      = h_so2        +h_x_ka1_so2    *chm1      +h_x_ka1_x_ka2_so2*chm1*chm1
       he_nh3      = h_nh3        +h_x_kb_kw_nh3  *ch
       he_hcooh    = h_hcooh      +h_x_ka_hcooh   *chm1
       he_ch3cooh  = h_ch3cooh    +h_x_ka_ch3cooh *chm1
       he_co2      = h_co2        +h_x_ka1_co2    *chm1      +h_x_ka1_x_ka2_co2*chm1*chm1
       
       !gas concentration in equilibrium with cloud water
       hno3g       = hno3    / (1. + he_hno3     * conv)
       so2g        = so2     / (1. + he_so2      * conv)
       nh3g        = nh3     / (1. + he_nh3      * conv)
       hcoohg      = hcooh   / (1. + he_hcooh    * conv)
       ch3coohg    = ch3cooh / (1. + he_ch3cooh  * conv)
       co2g        = co2     / (1. + he_co2      * conv)
       

       diagnostics(1) = h_x_ka_hno3    * hno3g    * patm / ch
       diagnostics(2) = h_x_ka_hcooh   * hcoohg   * patm / ch
       diagnostics(3) = h_x_ka_ch3cooh * ch3coohg  * patm / ch
       diagnostics(4) = h_x_ka1_so2    * so2g      *patm / ch
       diagnostics(5) = h_x_ka1_x_ka2_so2  * so2g  * patm / (ch*ch)
       diagnostics(6) = h_x_ka1_co2  * co2g        * patm / ch
       diagnostics(7) = h_x_ka1_x_ka2_co2  * co2g  * patm / (ch*ch)
       diagnostics(8) = e_so4
       diagnostics(9) = h_x_kb_kw_nh3  * nh3g     * patm * ch
       diagnostics(10) = kw/ch
       diagnostics(11) = sum(diagnostics(1:4)) + 2*diagnostics(5) + diagnostics(6) + 2*diagnostics(7)+2*diagnostics(8) - diagnostics(9) + diagnostics(10) - cH

    end if
        
  end subroutine cloud_pH_cubic

  subroutine cloud_pH_bisection(ihno3,iso2,inh3,ihcooh,ich3cooh,iso4,ico2,ialk,tk,patm,xlw,air_dens,ch,diagnostics)

    real*8, parameter      :: ch0 = 1.e-5       
    integer, parameter     :: iter_max = 100

    !input

    real*8, intent(in)     :: ihno3,iso2,inh3,ihcooh,ich3cooh ! vmr
    real*8, intent(in)     :: iso4,ico2, ialk                 ! vmr
    real*8, intent(in)     :: xlw                             ! l(water)/l(air)
    real*8, intent(in)     :: air_dens                        ! molec/cm3
    real*8, intent(in)     :: tk                              ! t[k]
    real*8, intent(in)     :: patm                            ! atmospheric pressure (atm)

    !ouput

    real*8, intent(out)    :: ch                              ! h concentration in mol/l
    real*8, intent(out),optional, dimension(cloud_nb_diag)  :: diagnostics

    !local variables
    real*8                 :: delta
    integer                :: iter
    logical                :: converged
    real*8                 :: kw
    real*8                 :: he_hno3,he_hcooh,he_ch3cooh,he_nh3, he_so2, he_co2
    real*8                 :: h_hno3,h_hcooh,h_ch3cooh,h_nh3, h_so2, h_co2
    real*8                 :: ka_hno3,ka_hcooh,ka_ch3cooh,kb_nh3, ka1_so2, ka2_so2, ka1_co2, ka2_co2
    real*8                 :: h_x_ka_hno3,h_x_ka1_so2,h_x_ka1_x_ka2_so2,h_x_kb_kw_nh3,h_x_ka_hcooh,h_x_ka_ch3cooh,h_x_ka1_co2,h_x_ka1_x_ka2_co2
    real*8                 :: e_hno3, e_so2_1,e_so2_2, e_so4, e_hcooh, e_ch3cooh, e_nh3, e_h2o, e_co2_1,e_co2_2,e_alk
    real*8                 :: hno3g, so2g, hcoohg, ch3coohg, nh3g, co2g
    real*8                 :: hno3, so2, hcooh, ch3cooh, nh3, co2, so4, alk
    real*8                 :: conv, chm1,com2

    !bisection
    real*8                 :: ph_min,ph_max,pH_current,neutral

    ch = ch0                      ! initial h concentration

    call so2aq(tk,h_so2,ka1_so2,ka2_so2)
    call hno3aq(tk,h_hno3,ka_hno3)
    call nh3aq(tk,h_nh3,kb_nh3)
    call hcoohaq(tk,h_hcooh,ka_hcooh)
    call ch3coohaq(tk,h_ch3cooh,ka_ch3cooh)
    call co2aq(tk,h_co2,ka1_co2,ka2_co2)
    call kw_aq(tk,kw)
    
    conv = ra*tk*xlw !-> atm/(mol/l(air))
    
    h_x_ka_hno3               = h_hno3*ka_hno3
    h_x_ka1_so2               = h_so2*ka1_so2
    h_x_ka1_x_ka2_so2         = h_so2*ka1_so2 * ka2_so2
    h_x_kb_kw_nh3             = h_nh3 * kb_nh3/kw
    h_x_ka_hcooh              = h_hcooh*ka_hcooh
    h_x_ka_ch3cooh            = h_ch3cooh * ka_ch3cooh
    h_x_ka1_co2               = h_co2*ka1_co2
    h_x_ka1_x_ka2_co2         = h_co2*ka1_co2*ka2_co2

    hno3    = max(ihno3,0.)
    so2     = max(iso2,0.)
    nh3     = max(inh3,0.)
    hcooh   = max(ihcooh,0.)
    ch3cooh = max(ich3cooh,0.)
    co2     = max(ico2,0.)
    so4     = max(iso4,0.)
    alk     = max(ialk,0.)
    
    !independent of ph
    e_so4 = 2.*so4*air_dens*const0/xlw ! vmr*molec/cm3*cm3/l*mole/molec * l(air)/l(water) -> vmr * mole/l(water)

    e_alk = alk*air_dens*const0/xlw

    iter  = 1

    ph_min = 1
    ph_max = 13
    
    do while ( iter .lt. iter_max )

       pH_current = (pH_min+ph_max)/2
       delta = cH
       cH = 10**(-pH_current)
       chm1 = 1./ch
       
       !calculate effective henry's constant          
       he_hno3     = h_hno3       +h_x_ka_hno3    *chm1
       he_so2      = h_so2        +h_x_ka1_so2    *chm1      +h_x_ka1_x_ka2_so2*chm1*chm1
       he_nh3      = h_nh3        +h_x_kb_kw_nh3  *ch
       he_hcooh    = h_hcooh      +h_x_ka_hcooh   *chm1
       he_ch3cooh  = h_ch3cooh    +h_x_ka_ch3cooh *chm1
       he_co2      = h_co2        +h_x_ka1_co2    *chm1      +h_x_ka1_x_ka2_co2*chm1*chm1
       
       !gas concentration in equilibrium with cloud water
       hno3g       = hno3    / (1. + he_hno3     * conv)
       so2g        = so2     / (1. + he_so2      * conv)
       nh3g        = nh3     / (1. + he_nh3      * conv)
       hcoohg      = hcooh   / (1. + he_hcooh    * conv)
       ch3coohg    = ch3cooh / (1. + he_ch3cooh  * conv)
       co2g        = co2     / (1. + he_co2      * conv)
       
       !concentrations in aqueous phase (equivalent)
       e_hno3     = h_x_ka_hno3    * hno3g    * patm
       e_so2_1    = h_x_ka1_so2    * so2g     * patm
       e_so2_2    = 2. * e_so2_1   * ka2_so2
       e_nh3      = h_x_kb_kw_nh3  * nh3g     * patm
       e_h2o      = kw
       e_co2_1    = h_x_ka1_co2     * co2g     * patm
       e_co2_2    = 2. * e_co2_1   * ka2_co2
       e_hcooh    = h_x_ka_hcooh   * hcoohg   * patm              
       e_ch3cooh  = h_x_ka_ch3cooh * ch3coohg * patm


       neutral  = e_hno3    * chm1 & 
                + e_hcooh   * chm1 & 
                + e_ch3cooh * chm1 &
                + e_so2_1   * chm1 &
                + e_so2_2   * chm1 * chm1  &
                + e_co2_1   * chm1         &
                + e_co2_2   * chm1 * chm1  &
                + e_so4                    &
                + kw        * chm1         &
                - e_nh3     * ch           &
                - e_alk                    &
                - ch

      delta = abs(ch-delta)/delta

       if ( abs( neutral / ch ) .gt. 1e-4 .and. delta .gt. 1e-4) then

          if (neutral <0) then
             pH_min = -log10(ch)
          else
             pH_max = -log10(ch)
          end if

          converged = .false.
       
       else
          converged = .true.
          exit
       end if
       
       iter = iter + 1
       
    end do



    if( .not. converged ) then
       write(*,*) 'pH not converged  (% change)=', &
            100.*abs(ch-delta)/delta
    end if


    if ( present(diagnostics) ) then

       chm1 = 1./ch
       
       !calculate effective henry's constant          
       he_hno3     = h_hno3       +h_x_ka_hno3    *chm1
       he_so2      = h_so2        +h_x_ka1_so2    *chm1      +h_x_ka1_x_ka2_so2*chm1*chm1
       he_nh3      = h_nh3        +h_x_kb_kw_nh3  *ch
       he_hcooh    = h_hcooh      +h_x_ka_hcooh   *chm1
       he_ch3cooh  = h_ch3cooh    +h_x_ka_ch3cooh *chm1
       he_co2      = h_co2        +h_x_ka1_co2    *chm1      +h_x_ka1_x_ka2_co2*chm1*chm1
       
       !gas concentration in equilibrium with cloud water
       hno3g       = hno3    / (1. + he_hno3     * conv)
       so2g        = so2     / (1. + he_so2      * conv)
       nh3g        = nh3     / (1. + he_nh3      * conv)
       hcoohg      = hcooh   / (1. + he_hcooh    * conv)
       ch3coohg    = ch3cooh / (1. + he_ch3cooh  * conv)
       co2g        = co2     / (1. + he_co2      * conv)
       
       diagnostics(1) = h_x_ka_hno3    * hno3g    * patm / ch
       diagnostics(2) = h_x_ka_hcooh   * hcoohg   * patm / ch
       diagnostics(3) = h_x_ka_ch3cooh * ch3coohg  * patm / ch
       diagnostics(4) = h_x_ka1_so2    * so2g      *patm / ch
       diagnostics(5) = h_x_ka1_x_ka2_so2  * so2g  * patm / (ch*ch)
       diagnostics(6) = h_x_ka1_co2  * co2g        * patm / ch
       diagnostics(7) = h_x_ka1_x_ka2_co2  * co2g  * patm / (ch*ch)
       diagnostics(8) = e_so4
       diagnostics(9) = h_x_kb_kw_nh3  * nh3g     * patm * ch
       diagnostics(10) = kw/ch
       diagnostics(11) = e_alk
       

    end if
        
  end subroutine cloud_pH_bisection



  subroutine henry_eff(n,tk,ch,h0,hv,tf,in_pa,he)

    real*8, intent(in)              :: ch, tk
    real*8, intent(in)              :: h0,hv,tf
    logical, intent(in)             :: in_pa
    integer                         :: n
    real*8, intent(out)             :: he
    logical                         :: error = .false.
    real*8                          :: ka,h,kb,kw,ka1,ka2

    if ( ch .lt. 0. ) then
       error = .true. 
    end if


    !calculate effective henry's constant          
    if ( n .eq. n_hno3 ) then
       call hno3aq(tk,h,ka)
       he     = h       * ( 1. + ka       *1./ch)
    else if (n .eq. n_so2 ) then
       call so2aq(tk,h,ka1,ka2)
       he      = h        *( 1. + ka1*1./ch      * (1. + ka2*1./ch) )
    else if (n .eq. n_nh3 ) then
       call nh3aq(tk,h,kb)
       call kw_aq(tk,kw)
       he      = h        *( 1. + kb/kw     *ch  )
    else if (n .eq. n_hcooh ) then
       call hcoohaq(tk,h,ka)
       he    = h      *( 1. + ka      *1./ch)
    else if (n .eq.n_ch3cooh ) then
       call ch3coohaq(tk,h,ka)
       he  = h    *( 1. + ka    *1./ch)
    else if (n .eq. n_h2o2 ) then
       call h2o2aq(tk,h,ka)
       he     = h       *( 1. + ka       *1./ch)
    else
       !use field table
       !this is in pa convert to atm
       he     = h0*exp(hv*tf)*1.01325e5
       error = .false.
    end if

    !all henry's constants are in mol/l/atm
    !some routines want h to be in mol/l/pa

    if (in_pa) then
       he = he/1.01325e5
    end if

    if (error) then
       write(*,*) 'Error: cH<0 in cloud_chem.F90'
    end if
    
    
  end subroutine henry_eff

  subroutine cloud_so2_chem(patm,ch,tk,xl,rso2_h2o2,rso2_o3,do_am3_bug)

    !all second order reactions have rates in l/mol//s

    real*8, intent(in)  :: ch
    real*8, intent(in)  :: tk
    real*8, intent(in)  :: xl !l(h2o)/l(air)
    real*8, intent(in)  :: patm
    logical, intent(in), optional  :: do_am3_bug

    real*8, intent( out) :: rso2_o3,rso2_h2o2
    

    !local variable
    real*8                :: conv,conv2
    real*8                :: h_so2,ka1_so2,ka2_so2
    real*8                :: h_o3,h_h2o2,ka_h2o2    
    real*8                :: he_so2,he_h2o2
    real*8                :: gfrac_so2,gfrac_h2o2,gfrac_o3
    real*8                :: kh2o2_hso3,ko3_so3,ko3_hso3,ko3_so2
    real*8                :: lfrac_so2,lfrac_hso3,lfrac_so3

    conv = ra * xl * tk ! atm/(mol/l)

    call so2aq(tk,h_so2,ka1_so2,ka2_so2)
    call o3aq(tk,h_o3)
    call h2o2aq(tk,h_h2o2,ka_h2o2)

    he_so2      = h_so2*(1+ka1_so2/ch*(1.+ka2_so2/ch))
    he_h2o2     = h_h2o2*(1+ka_h2o2/ch)

    !gas fraction
    gfrac_h2o2  = 1./(1.+he_h2o2 * conv)  !frac h2o2 in gas
    gfrac_so2   = 1./(1.+he_so2  * conv)  !frac so2  in gas
    gfrac_o3    = 1./(1.+h_o3    * conv)  !frac o3   in gas

!so2+h2o2

    !h2o2 + hso3
    kh2o2_hso3 =7.45e7*exp(-15.96*(298.15/tk-1.)) / (1.+13.*ch) !j05 

    rso2_h2o2  = kh2o2_hso3                           &
               * gfrac_h2o2 * h_h2o2                  &
               * gfrac_so2  * h_so2  * ka1_so2 
    
!so2+o3
    ko3_so2    = 2.4e4
    ko3_hso3   = 3.7e5 * exp( -18.56*(298.15/tk-1.))
    ko3_so3    = 1.5e9 * exp( -17.72*(298.15/tk-1.))

    lfrac_so2   = 1./(1.+ka1_so2/ch+ka2_so2*ka1_so2/(ch*ch))
    lfrac_hso3  = 1./(ch/ka1_so2 + 1. + ka2_so2/ch)
    lfrac_so3   = 1./(ch*ch/(ka1_so2*ka2_so2) + ch/ka2_so2 + 1.)

    rso2_o3    = gfrac_o3   * h_o3                &
               * gfrac_so2  * he_so2              &
               * (lfrac_so2*ko3_so2 + lfrac_hso3*ko3_hso3 + lfrac_so3*ko3_so3)


!at that stage, rates are in r_o3 and r_h2o2 are in M-1 s-1 * M/atm * M/atm = M/s 1/atm**2
!if we multiply by Patm^2 and then multiply by LWC * Ra * T / P ( l(w)/l(air)*atm/(mol/l(air)*K) * K / atm = l(w)/mol(air)
!result in vmr/s

! If do_am3_bug, then return rates in original units,
! and do (incorrect) units conversion in calling routine
if (.not. present(do_am3_bug)) then
   conv2 = patm*conv
   rso2_o3   = rso2_o3 * conv2
   rso2_h2o2 = rso2_h2o2 * conv2
else if (.not. do_am3_bug) then
   conv2 = patm*conv
   rso2_o3   = rso2_o3 * conv2
   rso2_h2o2 = rso2_h2o2 * conv2
end if
    
  end subroutine cloud_so2_chem

  subroutine so2aq(t,h,ka1,ka2)

    real*8, intent(in)  :: t
    real*8, intent(out) :: h,ka1,ka2

    real*8              :: tr

    tr = 298.15/t - 1.

    h     = 1.22     * exp(10.55*tr)             !j05
    ka1   = 1.71e-2  * exp(7.04*tr)              !j05
    ka2   = 5.99e-8  * exp(3.74*tr)              !j05

  end subroutine so2aq

  subroutine co2aq(t,h,ka1,ka2)

    real*8, intent(in)  :: t
    real*8, intent(out) :: h,ka1,ka2

    real*8              :: tr, ltr

    tr    = 298.15/t - 1.
    ltr   = -tr+log(298.15/t)

    h     = 3.41e-2*exp(8.19*tr-28.93*ltr)  !j05
    ka1   = 4.30e-7*exp(-30.8*tr+31.81*ltr) !j05
    ka2   = 4.68e-11*exp(-5.99*tr+38.84*ltr)
    
  end subroutine co2aq

  subroutine hcoohaq(t,h,ka)

    real*8, intent(in)  :: t
    real*8, intent(out) :: h,ka

    real*8              :: tr

    tr    = 298.15/t-1.

    h  = 5.39e3*exp(18.9*tr)   !j05
    ka = 1.86e-4*exp(-0.05*tr) !j05

  end subroutine hcoohaq
  
  subroutine ch3coohaq(t,h,ka)

    real*8, intent(in)  :: t
    real*8, intent(out) :: h,ka

    real*8              :: tr

    tr    = 298.15/t-1.

    h  = 8.6e3*exp(21.58*tr)  !j05
    ka = 1.75e-5*exp(0.1*tr)  !j05

  end subroutine ch3coohaq

  subroutine o3aq(t,h)

    real*8, intent(in)  :: t
    real*8, intent(out) :: h

    real*8              :: tr

    tr    = 298.15/t-1.

    h  = 1.13e-2*exp(7.72*tr)  !j05

  end subroutine o3aq

  subroutine h2o2aq(t,h,ka)

    real*8, intent(in)  :: t
    real*8, intent(out) :: h,ka

    real*8              :: tr

    tr    = 298.15/t-1.

    ka  = 2.2e-12*exp(-12.52*tr)  !j05
    h   = 7.45e4*exp(22.21*tr)    !j05

  end subroutine h2o2aq

  subroutine hno3aq(t,h,ka)

    real*8, intent(in)  :: t
    real*8, intent(out) :: h,ka

    real*8              :: tr
    real*8              :: ltr

    tr    = 298.15/t-1.
    ltr   = -tr+log(298.15/t)

    h    = 2.1e5
    ka   = 15.4  * exp(29.17 * tr - 16.83 * ltr)                    

  end subroutine hno3aq

  subroutine nh3aq(t,h,kb)

    real*8, intent(in)  :: t
    real*8, intent(out) :: h,kb

    real*8              :: tr,ltr

    tr    = 298.15/t-1.
    ltr   = -tr+log(298.15/t)

    h   = 5.76e1*exp(13.79*tr-5.39*ltr)    !j05
    kb  = 1.81e-5*exp(-1.5*tr+26.92*ltr)   !j05
    
  end subroutine nh3aq

  subroutine kw_aq(t,kw)

    real*8, intent(in)  :: t
    real*8, intent(out) :: kw
    real*8              :: tr,ltr

    tr    = 298.15/t-1.
    ltr   = -tr+log(298.15/t)

    kw=1.01e-14*exp(-22.52*tr+26.92*ltr)

  end subroutine kw_aq

  subroutine cubic( a2, a1, a0, nr, crutes )                                                                                                                                       
    integer           :: nr
    real*8              :: a2, a1, a0
    real*8              :: crutes(3)

! formulae can be found in numer. recip.  on page 145                         
!   kiran  developed  this version on 25/4/1990                                                                                                                                                                        

    real*8, parameter :: one    = 1.0d0
    real*8, parameter :: sqrt3  = 1.732050808d0
    real*8, parameter :: one3rd = 0.333333333d0
                    
    real*8            :: qq,    rr,    a2sq,  theta, dum1, dum2
    real*8            :: part1, part2, part3, rrsq,  phi,  yy1
    real*8            :: yy2,   yy3,   costh, sinth

    a2sq = a2 * a2
    qq   = ( a2sq - 3.d0*a1 ) / 9.d0
    rr   = ( a2*( 2.d0*a2sq - 9.d0*a1 ) + 27.d0*a0 ) / 54.d0

    dum1 = qq * qq * qq
    rrsq = rr * rr
    dum2 = dum1 - rrsq

    if ( dum2 .ge. 0.d0 ) then

       phi = sqrt( dum1 )

       if ( abs( phi ) .lt. 1.d-20 ) then
          crutes(1) = 0.0d0
          crutes(2) = 0.0d0
          crutes(3) = 0.0d0
          nr        = 0
       endif

       theta = acos( rr / phi ) / 3.0d0
       costh = cos( theta )
       sinth = sin( theta )
                 
       part1     = sqrt( qq )
       yy1       = part1 * costh
       yy2       = yy1 - a2/3.0d0
       yy3       = sqrt3 * part1 * sinth
       crutes(3) = -2.0d0*yy1 - a2/3.0d0
       crutes(2) = yy2 + yy3
       crutes(1) = yy2 - yy3
                 
       if ( crutes(1) .lt. 0.0d0 ) crutes(1) = 1.0d9
       if ( crutes(2) .lt. 0.0d0 ) crutes(2) = 1.0d9
       if ( crutes(3) .lt. 0.0d0 ) crutes(3) = 1.0d9
                 
       crutes(1) = min( crutes(1), crutes(2), crutes(3) )
       nr        = 3

    else

       part1     = sqrt( rrsq - dum1 )
       part2     = abs( rr )
       part3     = ( part1 + part2 )**one3rd
       crutes(1) = -sign(one,rr) * ( part3 + (qq/part3) ) - a2/3.d0
       crutes(2) = 0.d0
       crutes(3) = 0.d0
       nr        = 1

    endif

  end subroutine cubic


end module cloud_chem
