#ifdef CLUBB
! Alternate cloud fraction scheme

module alt_cloud_mod

use  sat_vapor_pres_mod,  only : compute_qs
use  constants_mod, only       : TFREEZE

implicit none
private

public :: alt_cloud

! Parameters
real, parameter :: rhmini       = 0.80          ! Minimum rh for ice cloud fraction > 0.
real, parameter :: rhmaxi       = 1.1           ! rhi at which ice cloud fraction = 1.
real, parameter :: qist_min     = 1.e-7         ! Minimum in-stratus ice IWC constraint [ kg/kg ]
real, parameter :: qist_max     = 5.e-3         ! Maximum in-stratus ice IWC constraint [ kg/kg ]
real, parameter :: minice       = 1.e-12   
real, parameter :: mincld       = 1.e-4   

contains

  subroutine alt_cloud(do_alt_cloud,dt,pfull,Tin,qvin,qain,qlin,qiin,tdt,qvdt,qadt,qldt,qidt)

  implicit none
  integer, intent(in)                   :: do_alt_cloud
  real, intent(in)                      :: dt
  real, dimension(:,:,:), intent(in)    :: pfull
  real, dimension(:,:,:), intent(in)    :: Tin,qvin,qain,qlin,qiin
  real, dimension(:,:,:), intent(inout) :: tdt,qvdt,qadt,qldt,qidt

  real, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) :: T,qv,qa,ql,qi,qc,qs

  real rhi, rhdif, aist, icimr

  integer i, j, k
  integer idim, jdim, kdim
  idim = size(T,1)
  jdim = size(T,2)
  kdim = size(T,3)

  ! Updated input fields

  T = Tin + tdt*dt
  qv = qvin + qvdt*dt
  qa = qain + qadt*dt
  ql = qlin + qldt*dt
  qi = qiin + qidt*dt

  ! Check for unphysical temperature values
  do k=1,kdim
   do j=1,jdim
    do i=1,idim

       if (T(i,j,k).lt.-150.0+273.15 .or. T(i,j,k).gt.90+273.15) then
         write(*,'(a,3i4,f6.2,2f12.6)') 'alt_cloud: bad temperature ',   &
          i,j,k,T(i,j,k),ql(i,j,k),qi(i,j,k)
       end if

    end do
   end do
  end do

  ! Compute saturation specific humidity
  call compute_qs( T, pfull, qs )

  ! Compute alternate cloud fraction
  do k=1,kdim
   do j=1,jdim
    do i=1,idim

      if (T(i,j,k).lt.TFREEZE) then

        ! set rh ice cloud fraction
        rhi = (qv(i,j,k)+qi(i,j,k))/qs(i,j,k)
        rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
        aist = min(1.0, max(rhdif,0.0)**2)

        ! limiter to remove empty cloud and ice with no cloud
        ! and set icecld fraction to mincld if ice exists
        if (qi(i,j,k).lt.minice) then
          aist=0.0
        else
          aist=max(mincld,aist)
        endif

        ! enforce limits on icimr
        if (qi(i,j,k).ge.minice) then
          icimr=qi(i,j,k)/aist

          !minimum
          if (icimr.lt.qist_min) then
            aist = max(0.0,min(1.0,qi(i,j,k)/qist_min))
          endif
          !maximum
          if (icimr.gt.qist_max) then
            aist = max(0.0,min(1.0,qi(i,j,k)/qist_max))
          endif

        ! cloud fraction should be largest of clubb, alternate ice
        qa(i,j,k) = max( aist, qa(i,j,k) )

        endif
      endif

    end do
   end do
  end do

  ! Update cloud fraction tendency
  do k=1,kdim
   do j=1,jdim
    do i=1,idim
      qadt(i,j,k) = (qa(i,j,k)-qain(i,j,k))/dt
    end do
   end do
  end do

  end subroutine alt_cloud

end module alt_cloud_mod

#endif
