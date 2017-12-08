      module MO_JPL_MOD

implicit none
character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
logical                       :: module_is_initialized = .false.

      CONTAINS

      subroutine JPL( rate, m, factor, ko, kinf,plnplv )
!-----------------------------------------------------------------
!        ... Calculate JPL troe rate
!-----------------------------------------------------------------
      
      implicit none

!-----------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------
      integer,intent(in)::    plnplv
      real, intent(in)  ::   factor
      real, intent(in)  ::   ko(plnplv)
      real, intent(in)  ::   kinf(plnplv)
      real, intent(in)  ::   m(plnplv)
      real, intent(out) ::   rate(plnplv)

!-----------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------
      real  ::     xpo( SIZE(rate) )
      

      xpo(:)  = ko(:) * m(:) / kinf(:)
      rate(:) = ko(:) / (1. + xpo(:))
      xpo(:)  = LOG10( xpo(:) )
      xpo(:)  = 1. / (1. + xpo(:)*xpo(:))
      rate(:) = rate(:) * factor**xpo(:)

      end subroutine JPL

      end module MO_JPL_MOD
