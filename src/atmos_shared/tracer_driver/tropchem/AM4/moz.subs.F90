





      module mo_setrxt_mod

      private
      public :: setrxt

      contains

      subroutine setrxt( rate, temp, m, plonl, plev, plnplv )


      use CHEM_MODS_MOD, only : rxntot
      use mo_jpl_mod, only : jpl

      implicit none

!-------------------------------------------------------
! ... Dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: plonl, plev, plnplv
      real, intent(in) :: temp(plonl,plev), m(plonl,plev)
      real, intent(inout) :: rate(plonl,plev,rxntot)

!-------------------------------------------------------
! ... Local variables
!-------------------------------------------------------
      real :: itemp(plonl,plev), exp_fac(plonl,plev)
      real, dimension(plonl,plev) :: ko, kinf

      rate(:,:,50) = 1.2e-10
      rate(:,:,56) = 1.8e-12
      rate(:,:,58) = 1.8e-12
      rate(:,:,69) = 3.5e-12
      rate(:,:,72) = 0.
      rate(:,:,80) = 1.5e-10
      rate(:,:,90) = 1.e-14
      rate(:,:,111) = 3.e-13
      rate(:,:,112) = 4.1e-14
      rate(:,:,128) = 8.37E-14
      rate(:,:,129) = 1.54E-13
      rate(:,:,133) = 3.7E-19
      rate(:,:,143) = 8.37E-14
      rate(:,:,145) = 1.6E-12
      rate(:,:,149) = 3.4E-15
      rate(:,:,152) = 8.37E-14
      rate(:,:,155) = 3.20E-12
      rate(:,:,167) = 8.37E-14
      rate(:,:,174) = 2.90E-11
      rate(:,:,175) = 1.00E-11
      rate(:,:,177) = 4.0E-16
      rate(:,:,178) = 1.50E-11
      rate(:,:,183) = 2.30E-12
      rate(:,:,188) = 4.00E-12
      rate(:,:,193) = 1.6e-12
      rate(:,:,247) = 3.17e-8
      itemp(:,:) = 1. / temp(:,:)
      rate(:,:,45) = 8e-12 * exp( -2060. * itemp(:,:) )
      rate(:,:,46) = 1.5e-11 * exp( -3600. * itemp(:,:) )
      rate(:,:,47) = 2.1e-11 * exp( 100. * itemp(:,:) )
      exp_fac(:,:) = exp( 180. * itemp(:,:) )
      rate(:,:,51) = 1.8e-11 * exp_fac(:,:)
      rate(:,:,89) = 4.2e-12 * exp_fac(:,:)
      rate(:,:,116) = 4.2e-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 200. * itemp(:,:) )
      rate(:,:,52) = 3e-11 * exp_fac(:,:)
      rate(:,:,85) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,98) = 8.78E-12 * exp_fac(:,:)
      rate(:,:,105) = 3.80e-12 * exp_fac(:,:)
      rate(:,:,113) = 5.18e-12 * exp_fac(:,:)
      rate(:,:,119) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,134) = 4.75E-12 * exp_fac(:,:)
      rate(:,:,144) = 8.78E-12 * exp_fac(:,:)
      rate(:,:,153) = 1.84E-12 * exp_fac(:,:)
      rate(:,:,164) = 6.13E-13 * exp_fac(:,:)
      rate(:,:,172) = 2.66e-12 * exp_fac(:,:)
      rate(:,:,173) = 1.14e-12 * exp_fac(:,:)
      rate(:,:,186) = 5.18E-12 * exp_fac(:,:)
      rate(:,:,231) = 5.5e-12 * exp_fac(:,:)
      rate(:,:,53) = 1.7e-12 * exp( -940. * itemp(:,:) )
      rate(:,:,54) = 1e-14 * exp( -490. * itemp(:,:) )
      rate(:,:,57) = 4.8e-11 * exp( 250. * itemp(:,:) )
      rate(:,:,59) = 2.8e-12 * exp( -1800. * itemp(:,:) )
      rate(:,:,60) = 2.15e-11 * exp( 110. * itemp(:,:) )
      rate(:,:,61) = 3.3e-11 * exp( 55. * itemp(:,:) )
      rate(:,:,62) = 1.63e-10 * exp( 60. * itemp(:,:) )
      exp_fac(:,:) = exp( 20. * itemp(:,:) )
      rate(:,:,63) = 7.25e-11 * exp_fac(:,:)
      rate(:,:,64) = 4.63e-11 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270. * itemp(:,:) )
      rate(:,:,65) = 3.3e-12 * exp_fac(:,:)
      rate(:,:,101) = 8.1e-12 * exp_fac(:,:)
      rate(:,:,218) = 7.4e-12 * exp_fac(:,:)
      rate(:,:,66) = 3e-12 * exp( -1500. * itemp(:,:) )
      rate(:,:,67) = 5.1e-12 * exp( 210. * itemp(:,:) )
      exp_fac(:,:) = exp( -2450. * itemp(:,:) )
      rate(:,:,68) = 1.2e-13 * exp_fac(:,:)
      rate(:,:,235) = 8.5e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 170. * itemp(:,:) )
      rate(:,:,75) = 1.5e-11 * exp_fac(:,:)
      rate(:,:,216) = 1.8e-11 * exp_fac(:,:)
      exp_fac(:,:) = exp( 380. * itemp(:,:) )
      rate(:,:,77) = 1.3e-12 * exp_fac(:,:)
      rate(:,:,130) = 3.61E-12 * exp_fac(:,:)
      rate(:,:,146) = 8.00E-12 * exp_fac(:,:)
      rate(:,:,154) = 4.40E-12 * exp_fac(:,:)
      rate(:,:,165) = 3.60E-12 * exp_fac(:,:)
      rate(:,:,79) = 2.45e-12 * exp( -1775. * itemp(:,:) )
      exp_fac(:,:) = exp( 300. * itemp(:,:) )
      rate(:,:,81) = 2.8e-12 * exp_fac(:,:)
      rate(:,:,170) = 2.80E-12 * exp_fac(:,:)
      rate(:,:,82) = 6.03e-13 * exp( -453. * itemp(:,:) )
      rate(:,:,83) = 2.30e-14 * exp( 677. * itemp(:,:) )
      rate(:,:,84) = 4.1e-13 * exp( 750. * itemp(:,:) )
      exp_fac(:,:) = exp( -1900. * itemp(:,:) )
      rate(:,:,86) = 3.4e-13 * exp_fac(:,:)
      rate(:,:,100) = 1.4e-12 * exp_fac(:,:)
      rate(:,:,87) = 5.5e-12 * exp( 125. * itemp(:,:) )
      rate(:,:,91) = 1.6e11 * exp( -4150. * itemp(:,:) )
      rate(:,:,92) = 1.2e-14 * exp( -2630. * itemp(:,:) )
      rate(:,:,94) = 5.50E-15 * exp( -1880. * itemp(:,:) )
      rate(:,:,95) = 4.6e-13 * exp( -1156. * itemp(:,:) )
      exp_fac(:,:) = exp( 350. * itemp(:,:) )
      rate(:,:,96) = 2.7e-12 * exp_fac(:,:)
      rate(:,:,99) = 4.63E-12 * exp_fac(:,:)
      rate(:,:,123) = 3.10E-11 * exp_fac(:,:)
      rate(:,:,126) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,138) = 2.70E-12 * exp_fac(:,:)
      rate(:,:,141) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,150) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,156) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,168) = 2.35E-12 * exp_fac(:,:)
      rate(:,:,169) = 0.35E-12 * exp_fac(:,:)
      rate(:,:,182) = 2.70E-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 700. * itemp(:,:) )
      rate(:,:,97) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,110) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,117) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,171) = 8.60E-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 980. * itemp(:,:) )
      rate(:,:,103) = 5.2e-13 * exp_fac(:,:)
      rate(:,:,190) = 5.20E-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500. * itemp(:,:) )
      rate(:,:,104) = 2.0e-12 * exp_fac(:,:)
      rate(:,:,107) = 2.5e-12 * exp_fac(:,:)
      rate(:,:,160) = 1.68E-12 * exp_fac(:,:)
      rate(:,:,161) = 1.87E-13 * exp_fac(:,:)
      rate(:,:,108) = 7.66e-12 * exp( -1020. * itemp(:,:) )
      rate(:,:,109) = 2.6e-12 * exp( 365. * itemp(:,:) )
      rate(:,:,114) = 1.55e-11 * exp( -540. * itemp(:,:) )
      rate(:,:,115) = 8.7e-12 * exp( -615. * itemp(:,:) )
      rate(:,:,118) = 3.75e-13 * exp( -40. * itemp(:,:) )
      rate(:,:,121) = 2.9e-12 * exp( -345. * itemp(:,:) )
      rate(:,:,122) = 6.9e-12 * exp( -230. * itemp(:,:) )
      rate(:,:,124) = 4.07E+08 * exp( -7694. * itemp(:,:) )
      rate(:,:,125) = 1.00e-14 * exp( -1970. * itemp(:,:) )
      exp_fac(:,:) = exp( 1300. * itemp(:,:) )
      rate(:,:,127) = 2.06E-13 * exp_fac(:,:)
      rate(:,:,137) = 2.06E-13 * exp_fac(:,:)
      rate(:,:,142) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,151) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,157) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,166) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,184) = 2.06E-13 * exp_fac(:,:)
      rate(:,:,131) = 2.4e-12 * exp( 360. * itemp(:,:) )
      rate(:,:,132) = 8.7e-14 * exp( 1650. * itemp(:,:) )
      exp_fac(:,:) = exp( 390. * itemp(:,:) )
      rate(:,:,135) = 1.90E-11 * exp_fac(:,:)
      rate(:,:,185) = 1.90E-11 * exp_fac(:,:)
      rate(:,:,136) = 5.78E-11 * exp( -400. * itemp(:,:) )
      rate(:,:,139) = 2.6E-12 * exp( 610. * itemp(:,:) )
      exp_fac(:,:) = exp( -1520. * itemp(:,:) )
      rate(:,:,140) = 8.5E-16 * exp_fac(:,:)
      rate(:,:,191) = 4.15E-15 * exp_fac(:,:)
      rate(:,:,147) = 2.90E+07 * exp( -5297. * itemp(:,:) )
      rate(:,:,148) = 1.40E-15 * exp( -2100. * itemp(:,:) )
      exp_fac(:,:) = exp( 340. * itemp(:,:) )
      rate(:,:,158) = 6.7E-12 * exp_fac(:,:)
      rate(:,:,176) = 3.1E-12 * exp_fac(:,:)
      rate(:,:,189) = 6.70E-12 * exp_fac(:,:)
      rate(:,:,159) = 4.3E-13 * exp( 1040. * itemp(:,:) )
      rate(:,:,179) = 1.40E-12 * exp( -1860. * itemp(:,:) )
      rate(:,:,180) = 1.60E-12 * exp( 305. * itemp(:,:) )
      rate(:,:,181) = 3.30E-12 * exp( -450. * itemp(:,:) )
      rate(:,:,187) = 3.15E-13 * exp( -448. * itemp(:,:) )
      rate(:,:,192) = 7.48E-12 * exp( 410. * itemp(:,:) )
      rate(:,:,194) = 1.2e-11 * exp( 440. * itemp(:,:) )
      rate(:,:,195) = 5.3e-16 * exp( -530. * itemp(:,:) )
      rate(:,:,196) = 1.2e-12 * exp( 490. * itemp(:,:) )
      rate(:,:,202) = 1.2e-11 * exp( -280. * itemp(:,:) )
      rate(:,:,204) = 1.90e-13 * exp( 530. * itemp(:,:) )
      rate(:,:,206) = 1.7e-12 * exp( -710. * itemp(:,:) )
      rate(:,:,207) = 1.4e-10 * exp( -470. * itemp(:,:) )
      rate(:,:,209) = 2.3e-11 * exp( -200. * itemp(:,:) )
      rate(:,:,210) = 2.8e-11 * exp( 85. * itemp(:,:) )
      exp_fac(:,:) = exp( 290. * itemp(:,:) )
      rate(:,:,211) = 6.4e-12 * exp_fac(:,:)
      rate(:,:,232) = 4.1e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( -800. * itemp(:,:) )
      rate(:,:,213) = 2.9e-12 * exp_fac(:,:)
      rate(:,:,223) = 1.7e-11 * exp_fac(:,:)
      rate(:,:,230) = 1.7e-11 * exp_fac(:,:)
      rate(:,:,214) = 7.3e-12 * exp( -1280. * itemp(:,:) )
      rate(:,:,215) = 2.6e-12 * exp( -350. * itemp(:,:) )
      exp_fac(:,:) = exp( 220. * itemp(:,:) )
      rate(:,:,217) = 2.7e-12 * exp_fac(:,:)
      rate(:,:,237) = 5.8e-12 * exp_fac(:,:)
      rate(:,:,219) = 8.1e-11 * exp( -30. * itemp(:,:) )
      exp_fac(:,:) = exp( 260. * itemp(:,:) )
      rate(:,:,225) = 2.3e-12 * exp_fac(:,:)
      rate(:,:,227) = 8.8e-12 * exp_fac(:,:)
      rate(:,:,226) = 4.5e-12 * exp( 460. * itemp(:,:) )
      rate(:,:,228) = 1.2e-10 * exp( -430. * itemp(:,:) )
      rate(:,:,229) = 4.8e-12 * exp( -310. * itemp(:,:) )
      rate(:,:,233) = 6.0e-13 * exp( 230. * itemp(:,:) )
      rate(:,:,234) = 4.5e-14 * exp( -1260. * itemp(:,:) )

      itemp(:,:) = 300. * itemp(:,:)

      ko(:,:) = 5.9e-33 * itemp(:,:)**1.4
      kinf(:,:) = 1.1e-12 * itemp(:,:)**-1.3
      call jpl( rate(1,1,48), m, .6, ko, kinf, plnplv )

      ko(:,:) = 1.5e-13 * itemp(:,:)**-0.6
      kinf(:,:) = 2.1e9 * itemp(:,:)**-6.1
      call jpl( rate(1,1,49), m, .6, ko, kinf, plnplv )

      ko(:,:) = 2.e-30 * itemp(:,:)**4.4
      kinf(:,:) = 1.4e-12 * itemp(:,:)**.7
      call jpl( rate(1,1,70), m, .6, ko, kinf, plnplv )

      ko(:,:) = 1.8e-30 * itemp(:,:)**3.0
      kinf(:,:) = 2.8e-11
      call jpl( rate(1,1,73), m, .6, ko, kinf, plnplv )

      ko(:,:) = 2.0e-31 * itemp(:,:)**3.4
      kinf(:,:) = 2.9e-12 * itemp(:,:)**1.1
      call jpl( rate(1,1,76), m, .6, ko, kinf, plnplv )

      ko(:,:) = 1.e-28 * itemp(:,:)**4.5
      kinf(:,:) = 7.5e-12 * itemp(:,:)**0.85
      call jpl( rate(1,1,88), m, .6, ko, kinf, plnplv )

      ko(:,:) = 8.e-27 * itemp(:,:)**3.5
      kinf(:,:) = 3.e-11
      call jpl( rate(1,1,93), m, .5, ko, kinf, plnplv )

      ko(:,:) = 9.7e-29 * itemp(:,:)**5.6
      kinf(:,:) = 9.3e-12 * itemp(:,:)**1.5
      call jpl( rate(1,1,102), m, .6, ko, kinf, plnplv )

      ko(:,:) = 9.0e-28 * itemp(:,:)**8.9
      kinf(:,:) = 7.7e-12 * itemp(:,:)**.2
      call jpl( rate(1,1,162), m, .6, ko, kinf, plnplv )

      ko(:,:) = 3.3e-31 * itemp(:,:)**4.3
      kinf(:,:) = 1.6e-12
      call jpl( rate(1,1,201), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 4.4e-32 * itemp(:,:)**1.3
      kinf(:,:) = 4.7e-11 * itemp(:,:)**0.2
      call jpl( rate(1,1,208), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 1.8e-31 * itemp(:,:)**3.4
      kinf(:,:) = 1.5e-11 * itemp(:,:)**1.9
      call jpl( rate(1,1,212), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 6.9e-31 * itemp(:,:)**1.0
      kinf(:,:) = 2.6e-11
      call jpl( rate(1,1,220), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 1.6e-32 * itemp(:,:)**4.5
      kinf(:,:) = 2.0e-12 * itemp(:,:)**2.4
      call jpl( rate(1,1,221), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 5.2e-31 * itemp(:,:)**3.2
      kinf(:,:) = 6.9e-12 * itemp(:,:)**2.9
      call jpl( rate(1,1,224), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 9.0e-32 * itemp(:,:)**1.5
      kinf(:,:) = 3.0e-11
      call jpl( rate(1,1,236), m, 0.6, ko, kinf, plnplv )

      end subroutine setrxt

      end module mo_setrxt_mod

      module mo_adjrxt_mod

      private
      public :: adjrxt

      contains

      subroutine adjrxt( rate, inv, m, plnplv )

      use CHEM_MODS_MOD, only : nfs, rxntot

      implicit none

!--------------------------------------------------------------------
! ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(in) :: inv(plnplv,nfs)
      real, intent(in) :: m(plnplv)
      real, intent(inout) :: rate(plnplv,rxntot)

!--------------------------------------------------------------------
! ... Local variables
!--------------------------------------------------------------------
      real :: im(plnplv)

      rate(:, 46) = rate(:, 46) * inv(:, 3)
      rate(:, 48) = rate(:, 48) * inv(:, 1)
      rate(:, 60) = rate(:, 60) * inv(:, 2)
      rate(:, 61) = rate(:, 61) * inv(:, 3)
      rate(:, 70) = rate(:, 70) * inv(:, 1)
      rate(:, 71) = rate(:, 71) * inv(:, 1)
      rate(:, 73) = rate(:, 73) * inv(:, 1)
      rate(:, 76) = rate(:, 76) * inv(:, 1)
      rate(:, 78) = rate(:, 78) * inv(:, 1)
      rate(:, 88) = rate(:, 88) * inv(:, 1)
      rate(:, 90) = rate(:, 90) * inv(:, 3)
      rate(:, 93) = rate(:, 93) * inv(:, 1)
      rate(:,102) = rate(:,102) * inv(:, 1)
      rate(:,106) = rate(:,106) * inv(:, 1)
      rate(:,162) = rate(:,162) * inv(:, 1)
      rate(:,201) = rate(:,201) * inv(:, 1)
      rate(:,212) = rate(:,212) * inv(:, 1)
      rate(:,220) = rate(:,220) * inv(:, 1)
      rate(:,221) = rate(:,221) * inv(:, 1)
      rate(:,222) = rate(:,222) * inv(:, 1)
      rate(:,224) = rate(:,224) * inv(:, 1)
      rate(:,236) = rate(:,236) * inv(:, 1)
      rate(:, 44) = rate(:, 44) * inv(:, 3) * inv(:, 1)
      rate(:,208) = rate(:,208) * inv(:, 3) * inv(:, 1)
      rate(:, 45) = rate(:, 45) * m(:)
      rate(:, 47) = rate(:, 47) * m(:)
      rate(:, 48) = rate(:, 48) * m(:)
      rate(:, 49) = rate(:, 49) * m(:)
      rate(:, 50) = rate(:, 50) * m(:)
      rate(:, 51) = rate(:, 51) * m(:)
      rate(:, 52) = rate(:, 52) * m(:)
      rate(:, 53) = rate(:, 53) * m(:)
      rate(:, 54) = rate(:, 54) * m(:)
      rate(:, 55) = rate(:, 55) * m(:)
      rate(:, 56) = rate(:, 56) * m(:)
      rate(:, 57) = rate(:, 57) * m(:)
      rate(:, 58) = rate(:, 58) * m(:)
      rate(:, 59) = rate(:, 59) * m(:)
      rate(:, 62) = rate(:, 62) * m(:)
      rate(:, 63) = rate(:, 63) * m(:)
      rate(:, 64) = rate(:, 64) * m(:)
      rate(:, 65) = rate(:, 65) * m(:)
      rate(:, 66) = rate(:, 66) * m(:)
      rate(:, 67) = rate(:, 67) * m(:)
      rate(:, 68) = rate(:, 68) * m(:)
      rate(:, 69) = rate(:, 69) * m(:)
      rate(:, 70) = rate(:, 70) * m(:)
      rate(:, 72) = rate(:, 72) * m(:)
      rate(:, 73) = rate(:, 73) * m(:)
      rate(:, 74) = rate(:, 74) * m(:)
      rate(:, 75) = rate(:, 75) * m(:)
      rate(:, 76) = rate(:, 76) * m(:)
      rate(:, 77) = rate(:, 77) * m(:)
      rate(:, 79) = rate(:, 79) * m(:)
      rate(:, 80) = rate(:, 80) * m(:)
      rate(:, 81) = rate(:, 81) * m(:)
      rate(:, 82) = rate(:, 82) * m(:)
      rate(:, 83) = rate(:, 83) * m(:)
      rate(:, 84) = rate(:, 84) * m(:)
      rate(:, 85) = rate(:, 85) * m(:)
      rate(:, 86) = rate(:, 86) * m(:)
      rate(:, 87) = rate(:, 87) * m(:)
      rate(:, 88) = rate(:, 88) * m(:)
      rate(:, 89) = rate(:, 89) * m(:)
      rate(:, 92) = rate(:, 92) * m(:)
      rate(:, 93) = rate(:, 93) * m(:)
      rate(:, 94) = rate(:, 94) * m(:)
      rate(:, 95) = rate(:, 95) * m(:)
      rate(:, 96) = rate(:, 96) * m(:)
      rate(:, 97) = rate(:, 97) * m(:)
      rate(:, 98) = rate(:, 98) * m(:)
      rate(:, 99) = rate(:, 99) * m(:)
      rate(:,100) = rate(:,100) * m(:)
      rate(:,101) = rate(:,101) * m(:)
      rate(:,102) = rate(:,102) * m(:)
      rate(:,103) = rate(:,103) * m(:)
      rate(:,104) = rate(:,104) * m(:)
      rate(:,105) = rate(:,105) * m(:)
      rate(:,107) = rate(:,107) * m(:)
      rate(:,108) = rate(:,108) * m(:)
      rate(:,109) = rate(:,109) * m(:)
      rate(:,110) = rate(:,110) * m(:)
      rate(:,111) = rate(:,111) * m(:)
      rate(:,112) = rate(:,112) * m(:)
      rate(:,113) = rate(:,113) * m(:)
      rate(:,114) = rate(:,114) * m(:)
      rate(:,115) = rate(:,115) * m(:)
      rate(:,116) = rate(:,116) * m(:)
      rate(:,117) = rate(:,117) * m(:)
      rate(:,118) = rate(:,118) * m(:)
      rate(:,119) = rate(:,119) * m(:)
      rate(:,120) = rate(:,120) * m(:)
      rate(:,121) = rate(:,121) * m(:)
      rate(:,122) = rate(:,122) * m(:)
      rate(:,123) = rate(:,123) * m(:)
      rate(:,125) = rate(:,125) * m(:)
      rate(:,126) = rate(:,126) * m(:)
      rate(:,127) = rate(:,127) * m(:)
      rate(:,128) = rate(:,128) * m(:)
      rate(:,129) = rate(:,129) * m(:)
      rate(:,130) = rate(:,130) * m(:)
      rate(:,131) = rate(:,131) * m(:)
      rate(:,132) = rate(:,132) * m(:)
      rate(:,133) = rate(:,133) * m(:)
      rate(:,134) = rate(:,134) * m(:)
      rate(:,135) = rate(:,135) * m(:)
      rate(:,136) = rate(:,136) * m(:)
      rate(:,137) = rate(:,137) * m(:)
      rate(:,138) = rate(:,138) * m(:)
      rate(:,139) = rate(:,139) * m(:)
      rate(:,140) = rate(:,140) * m(:)
      rate(:,141) = rate(:,141) * m(:)
      rate(:,142) = rate(:,142) * m(:)
      rate(:,143) = rate(:,143) * m(:)
      rate(:,144) = rate(:,144) * m(:)
      rate(:,145) = rate(:,145) * m(:)
      rate(:,146) = rate(:,146) * m(:)
      rate(:,148) = rate(:,148) * m(:)
      rate(:,149) = rate(:,149) * m(:)
      rate(:,150) = rate(:,150) * m(:)
      rate(:,151) = rate(:,151) * m(:)
      rate(:,152) = rate(:,152) * m(:)
      rate(:,153) = rate(:,153) * m(:)
      rate(:,154) = rate(:,154) * m(:)
      rate(:,155) = rate(:,155) * m(:)
      rate(:,156) = rate(:,156) * m(:)
      rate(:,157) = rate(:,157) * m(:)
      rate(:,158) = rate(:,158) * m(:)
      rate(:,159) = rate(:,159) * m(:)
      rate(:,160) = rate(:,160) * m(:)
      rate(:,161) = rate(:,161) * m(:)
      rate(:,162) = rate(:,162) * m(:)
      rate(:,164) = rate(:,164) * m(:)
      rate(:,165) = rate(:,165) * m(:)
      rate(:,166) = rate(:,166) * m(:)
      rate(:,167) = rate(:,167) * m(:)
      rate(:,168) = rate(:,168) * m(:)
      rate(:,169) = rate(:,169) * m(:)
      rate(:,170) = rate(:,170) * m(:)
      rate(:,171) = rate(:,171) * m(:)
      rate(:,172) = rate(:,172) * m(:)
      rate(:,173) = rate(:,173) * m(:)
      rate(:,174) = rate(:,174) * m(:)
      rate(:,175) = rate(:,175) * m(:)
      rate(:,176) = rate(:,176) * m(:)
      rate(:,177) = rate(:,177) * m(:)
      rate(:,178) = rate(:,178) * m(:)
      rate(:,179) = rate(:,179) * m(:)
      rate(:,180) = rate(:,180) * m(:)
      rate(:,181) = rate(:,181) * m(:)
      rate(:,182) = rate(:,182) * m(:)
      rate(:,183) = rate(:,183) * m(:)
      rate(:,184) = rate(:,184) * m(:)
      rate(:,185) = rate(:,185) * m(:)
      rate(:,186) = rate(:,186) * m(:)
      rate(:,187) = rate(:,187) * m(:)
      rate(:,188) = rate(:,188) * m(:)
      rate(:,189) = rate(:,189) * m(:)
      rate(:,190) = rate(:,190) * m(:)
      rate(:,191) = rate(:,191) * m(:)
      rate(:,192) = rate(:,192) * m(:)
      rate(:,193) = rate(:,193) * m(:)
      rate(:,194) = rate(:,194) * m(:)
      rate(:,195) = rate(:,195) * m(:)
      rate(:,196) = rate(:,196) * m(:)
      rate(:,201) = rate(:,201) * m(:)
      rate(:,202) = rate(:,202) * m(:)
      rate(:,203) = rate(:,203) * m(:)
      rate(:,204) = rate(:,204) * m(:)
      rate(:,206) = rate(:,206) * m(:)
      rate(:,207) = rate(:,207) * m(:)
      rate(:,209) = rate(:,209) * m(:)
      rate(:,210) = rate(:,210) * m(:)
      rate(:,211) = rate(:,211) * m(:)
      rate(:,212) = rate(:,212) * m(:)
      rate(:,213) = rate(:,213) * m(:)
      rate(:,214) = rate(:,214) * m(:)
      rate(:,215) = rate(:,215) * m(:)
      rate(:,216) = rate(:,216) * m(:)
      rate(:,217) = rate(:,217) * m(:)
      rate(:,218) = rate(:,218) * m(:)
      rate(:,219) = rate(:,219) * m(:)
      rate(:,220) = rate(:,220) * m(:)
      rate(:,221) = rate(:,221) * m(:)
      rate(:,223) = rate(:,223) * m(:)
      rate(:,224) = rate(:,224) * m(:)
      rate(:,225) = rate(:,225) * m(:)
      rate(:,226) = rate(:,226) * m(:)
      rate(:,227) = rate(:,227) * m(:)
      rate(:,228) = rate(:,228) * m(:)
      rate(:,229) = rate(:,229) * m(:)
      rate(:,230) = rate(:,230) * m(:)
      rate(:,231) = rate(:,231) * m(:)
      rate(:,232) = rate(:,232) * m(:)
      rate(:,233) = rate(:,233) * m(:)
      rate(:,234) = rate(:,234) * m(:)
      rate(:,235) = rate(:,235) * m(:)
      rate(:,236) = rate(:,236) * m(:)
      rate(:,237) = rate(:,237) * m(:)
      rate(:,238) = rate(:,238) * m(:)
      rate(:,239) = rate(:,239) * m(:)
      rate(:,240) = rate(:,240) * m(:)
      rate(:,241) = rate(:,241) * m(:)
      rate(:,242) = rate(:,242) * m(:)
      rate(:,243) = rate(:,243) * m(:)
      rate(:,244) = rate(:,244) * m(:)
      rate(:,245) = rate(:,245) * m(:)
      rate(:,246) = rate(:,246) * m(:)

      end subroutine adjrxt

      end module mo_adjrxt_mod

      module mo_phtadj_mod

      private
      public :: phtadj

      contains

      subroutine phtadj( p_rate, inv, m, plnplv )

      use CHEM_MODS_MOD, only : nfs, phtcnt

      implicit none

!--------------------------------------------------------------------
! ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(in) :: inv(plnplv,nfs)
      real, intent(in) :: m(plnplv)
      real, intent(inout) :: p_rate(plnplv,phtcnt)

!--------------------------------------------------------------------
! ... Local variables
!--------------------------------------------------------------------
      real :: im(plnplv)

      im(:) = 1. / m(:)
      p_rate(:, 1) = p_rate(:, 1) * inv(:, 3) * im(:)

      end subroutine phtadj

      end module mo_phtadj_mod

      module mo_rxt_mod

      private
      public :: rxt_mod

      contains

      subroutine rxt_mod( rate, het_rates, grp_ratios, plnplv )

      use CHEM_MODS_MOD, only : rxntot, hetcnt, grpcnt

      implicit none

!---------------------------------------------------------------------------
! ... Dummy arguments
!---------------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(inout) :: rate(plnplv,rxntot)
      real, intent(inout) :: het_rates(plnplv,hetcnt)
      real, intent(in) :: grp_ratios(plnplv,grpcnt)


      end subroutine rxt_mod

      end module mo_rxt_mod

      module mo_make_grp_vmr_mod

      private
      public :: mak_grp_vmr

      contains

      subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, plonl )

      use MO_GRID_MOD, only : plev, pcnstm1
      use CHEM_MODS_MOD, only : grpcnt

      implicit none

!----------------------------------------------------------------------------
! ... Dummy arguments
!----------------------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in) :: vmr(plonl,plev,pcnstm1)
      real, intent(in) :: group_ratios(plonl,plev,grpcnt)
      real, intent(out) :: group_vmrs(plonl,plev,grpcnt)

!----------------------------------------------------------------------------
! ... Local variables
!----------------------------------------------------------------------------
      integer :: k

      end subroutine mak_grp_vmr

      end module mo_make_grp_vmr_mod
