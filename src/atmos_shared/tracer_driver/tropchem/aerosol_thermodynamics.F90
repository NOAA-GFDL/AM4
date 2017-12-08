module aerosol_thermodynamics

  use tropchem_types_mod, only : small_value

  private
  public :: aerosol_thermo
  public :: AERO_ISORROPIA, AERO_LEGACY, NO_AERO

  integer :: AERO_LEGACY    = 1
  integer :: AERO_ISORROPIA = 2
  integer :: NO_AERO        = 3

contains

  !do_aerosol_thermo
  subroutine aerosol_thermo( aerosol_thermo_type, rh, TK, pres, so4, nh3, nh4, hno3, no3)

    real, intent(in) :: rh   ! relative humidity (unitles)
    real, intent(in) :: TK   ! temperature (K)
    real, intent(in) :: pres ! pressure (Pa)
    real, intent(in) :: so4  ! sulfate vmr
    real, intent(inout) :: nh3  ! gas-phase ammonia vmr
    real, intent(inout) :: nh4  ! ammonium sulfate/bisulfate vmr
    real, intent(inout) :: no3  ! ammonium nitrate vmr
    real, intent(inout) :: hno3 ! nitric acid vmr
    integer, intent(in) :: aerosol_thermo_type

    if (aerosol_thermo_type .eq. AERO_ISORROPIA) then
       call isorropia_interface(rh, TK, pres, so4, nh3, nh4, hno3, no3)
    elseif ( aerosol_thermo_type .eq. AERO_LEGACY ) then
       call old_am3_interface(rh,TK,nh3,hno3,no3)
    end if

  end subroutine aerosol_thermo

  !do_old_am3
  subroutine old_am3_interface(rh,tz,xnh3,xhno3,xant)

    real, intent(in)  :: rh
    real, intent(inout)  :: xnh3, xhno3, xant
    real, intent(in) :: tz !T in K
    real              :: xx0, yy1, com, com1, xkp, xra
    real              :: cnh3, chno3

    real, parameter ::  xa0 = 11.,   &
         xb0 = -.1,   &
         xa1 = 1.053, &
         xb1 = -4.368,&
         xa2 = 1.016, &
         xb2 = -2.54, &
         xa3 = .816e-32, &
         xb3 = .259

    xx0 = xa0 + xb0*RH

    if( RH >= 90. ) then
       yy1 = xa1*EXP( xb1/xx0 )
    else
       yy1 = xa2*EXP( xb2/xx0 )
    end if

    xkp = yy1*(xa3*EXP( xb3*tz )/.7) &    ! ppb**2
         * 1.e-18                      ! mixing ratio

    cnh3  = xnh3
    chno3 = xhno3
    com = cnh3*chno3

    com1 = (cnh3 + chno3)**2 - 4.*(cnh3*chno3 - xkp)
    com1 = MAX( com1,1.e-30 )

    if( com >= xkp ) then   ! NH4NO3 is formed
       xra       = .5*(cnh3 + chno3 - SQRT(com1))
       xant      = MAX( xant + xra, small_value )
       xnh3      = MAX( xnh3 - xra, small_value )
       xhno3     = MAX( xhno3- xra, small_value )
    end if

  end subroutine old_am3_interface

  !do_isorropia
  subroutine isorropia_interface( rh, temp, pres, so4, nh3, nh4, hno3, no3)

    real, intent(in)       :: rh   ! relative humidity (unitles)
    real, intent(in)       :: pres ! pressure in Pa
    real, intent(in)       :: so4  ! sulfate vmr
    real, intent(inout)    :: nh3  ! gas-phase ammonia vmr
    real, intent(inout)    :: nh4  ! ammonium sulfate vmr
    real, intent(inout)    :: no3  ! ammonium nitrate vmr
    real, intent(inout)    :: hno3 ! nitric acid
    real, intent(in)       :: temp ! K


    integer, parameter     :: nothera =  9
    integer, parameter     :: nctrla  =  2
    integer, parameter     :: ncompa  =  8
    integer, parameter     :: nionsa  = 10
    integer, parameter     :: ngasaqa =  3
    integer, parameter     :: nsldsa  = 19

    real                   :: tso4        ! mol/m3
    real                   :: tnh3        ! mol/m3
    real                   :: tno3        ! mol/m3
    real                   :: tna,tcl,tk,tmg,tca        ! mol/m3
    real                   :: air_den     ! mol/m3
    real                   :: inv_air_den ! m3/mol

    real                 :: aerliq(nionsa+ngasaqa+2)
    real                 :: aersld(nsldsa) 
    real                 :: gas(ngasaqa) 
    real                 :: other(nothera)
    real                 :: wi(ncompa)    
    real                 :: wt(ncompa)
    real                 :: cntrl(nctrla)
    real                 :: rhi
    real                 :: ratio_nhx,ratio_no3
    character(len=15)    :: scasi

    !fwd problem only 
    cntrl(1) = 0.0 

    ! Metastable for now 
    cntrl(2) = 1.0

    !air density in mol/m3
    air_den = pres/(8.314e0*temp)
    inv_air_den = 1e0/air_den

    tso4 = so4*air_den
    tno3 = (no3+hno3)*air_den
    tnh3 = (nh3+nh4)*air_den

    !for now no na, cl, mg, ca, k
    tna  = 0.e0
    tcl  = 0.e0
    tca  = 0.e0
    tk   = 0.e0
    tmg  = 0.e0

    !setup initial arrays for isorropia
    !insert concentrations [mole/m3] into wi & prevent underflow        
    wi(1)    = max( tna,  small_value )
    wi(2)    = max( tso4, small_value )
    wi(3)    = max( tnh3, small_value )
    wi(4)    = max( tno3, small_value )
    wi(5)    = max( tcl,  small_value )
    wi(6)    = max( tca,  small_value )
    wi(7)    = max( tk,   small_value )
    wi(8)    = max( tmg,  small_value )

    !relative humidity
    rhi      = rh/100
    rhi      = max( 0.01,  rhi )
    rhi      = min( 0.995, rhi )

    ! perform aerosol thermodynamic equilibrium 
    ! isoropia can be found in isoropiaiicode.f
    ! inputs are wi, rhi, tempi, cntrl
    call isoropia (wi, rhi, temp, cntrl,        &
         wt, gas, aerliq, aersld,     &
         scasi, other)

!    if ( ltracer_nit ) then

       nh3      = max(min(gas(1),tnh3),0.)
       hno3     = max(min(gas(2),tno3),0.)
       no3      = max(tno3-hno3,0.)
       nh4      = max(tnh3-nh3,0.)

    ! else

    !    !here we need to be careful
    !    !in the stratosphere hno3 is predicted to go in the particle but there is no ammonia, this creates a net source of ammonia

    !    nh3      = max(min(gas(1),tnh3),0.)
    !    !if nitrate > (tnh3 -2*so4)
    !    no3      = max(min(tno3-gas(2), tnh3-2*tso4-nh3),0.) !no3 cannot exceed nh4-2*so4
    !    hno3 = max(tno3-no3,0.)
    !    nh4  = max(tnh3-nh3-no3,0.)

    ! end if


    hno3     = max( hno3 * inv_air_den , 0. )
    no3      = max( no3  * inv_air_den , 0. )
    nh4      = max( nh4  * inv_air_den,  0. )
    nh3      = max( nh3  * inv_air_den , 0. )

  end subroutine isorropia_interface

end module aerosol_thermodynamics
