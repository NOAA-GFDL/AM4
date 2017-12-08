module physics_radiation_exch_mod
#include <fms_platform.h>
!---------------------------------------------------------------------

use fms_mod,            only: lowercase, error_mesg, FATAL

!---- module data ----
use block_control_mod,  only: block_control_type

!---- public data ----
!---
!---Exch_ctrl
 public exchange_control_type
 type  exchange_control_type
     logical           :: doing_prog_clouds
     logical           :: doing_donner
     logical           :: doing_uw_conv
     logical           :: do_cosp
     logical           :: do_modis_yim
     logical           :: donner_meso_is_largescale
     integer           :: ncol
     integer           :: ncld
     real              :: min_diam_ice
     real              :: dcs
     real              :: min_diam_drop
     real              :: max_diam_drop
     real              :: cosp_frequency
     character(len=16) :: cloud_type_form  ! indicator of radiatively active clouds
     character(len=16) :: cosp_precip_sources ! indicator of precipitation sources that are seen in COSP
     real              :: qmin
     real              :: N_min
     real              :: N_land
     real              :: N_ocean
     real              :: qcvar
     integer           :: overlap
     logical           :: do_liq_num
     logical           :: do_ice_num
     integer           :: do_clubb
 end type  exchange_control_type

!---
!--- Radiation Flux block type
 public radiation_flux_block_type
 type radiation_flux_block_type
     real, dimension(:,:,:), _ALLOCATABLE :: tdt_rad           _NULL, &
                                             tdt_lw            _NULL, &
                                             extinction        _NULL
     real, dimension(:,:),   _ALLOCATABLE :: flux_sw           _NULL, &
                                             flux_sw_dir            _NULL, &
                                             flux_sw_dif            _NULL, &
                                             flux_sw_down_vis_dir   _NULL, &
                                             flux_sw_down_vis_dif   _NULL, &
                                             flux_sw_down_total_dir _NULL, &
                                             flux_sw_down_total_dif _NULL, &
                                             flux_sw_vis            _NULL, &
                                             flux_sw_vis_dir        _NULL, & 
                                             flux_sw_vis_dif        _NULL, &
                                             flux_lw                _NULL, &
                                             coszen                 _NULL
     contains
         procedure :: alloc=>alloc_radiation_flux_block_type
         procedure :: dealloc=>dealloc_radiation_flux_block_type
 end type radiation_flux_block_type

!--- Radiation Flux control type
 public radiation_flux_control_type
 type radiation_flux_control_type
     logical :: do_rad
 end type radiation_flux_control_type

!--- Rad_flux type definition
 public radiation_flux_type
 type radiation_flux_type
     type (radiation_flux_control_type) :: control
     type (radiation_flux_block_type), dimension(:), _ALLOCATABLE :: block _NULL
 end type radiation_flux_type


!---
!--- Cloud Scheme data type
public cloud_scheme_data_type
type cloud_scheme_data_type
     character(len=16) :: scheme_name
     real,  dimension(:,:,:), _ALLOCATABLE :: cloud_area     _NULL, & ! cell, meso, lsc, shallow
                                              liquid_amt     _NULL, & ! cell, meso, lsc, shallow
                                              ice_amt        _NULL, & ! cell, meso, lsc, shallow
                                              droplet_number _NULL, & ! cell, meso, lsc, shallow
                                              ice_number     _NULL, & ! lsc, shallow
                                              rain           _NULL, & ! lsc
                                              snow           _NULL, & ! lsc
                                              rain_size      _NULL, & ! lsc
                                              snow_size      _NULL, & ! lsc
                                              liquid_size    _NULL, & ! cell, meso
                                              ice_size       _NULL    ! cell, meso
     integer, dimension(:,:), _ALLOCATABLE :: nsum_out       _NULL    ! cell, meso
end type cloud_scheme_data_type

!--- Moist Clouds block type
 public clouds_from_moist_block_type
 type clouds_from_moist_block_type
     integer           :: index_strat
     integer           :: index_donner_meso
     integer           :: index_donner_cell
     integer           :: index_uw_conv
     type(cloud_scheme_data_type), dimension(:), _ALLOCATABLE :: Cloud_data _NULL
 end type clouds_from_moist_block_type

!--- Moist_clouds type definition
 public clouds_from_moist_type
 type clouds_from_moist_type
    type (clouds_from_moist_block_type), dimension(:), _ALLOCATABLE :: block _NULL
 end type clouds_from_moist_type


!---
!--- Cosp Rad block type
 public cosp_from_rad_block_type
 type cosp_from_rad_block_type
     real, dimension(:,:,:,:), _ALLOCATABLE :: tau_stoch        _NULL, &
                                               lwem_stoch       _NULL, &
                                               stoch_cloud_type _NULL, &
                                               stoch_conc_drop  _NULL, &
                                               stoch_conc_ice   _NULL, &
                                               stoch_size_drop  _NULL, &
                                               stoch_size_ice   _NULL
     real, dimension(:,:,:),   _ALLOCATABLE :: mr_ozone         _NULL
     real, dimension(:,:),     _ALLOCATABLE :: daytime          _NULL
     real, dimension(:,:),     _ALLOCATABLE :: tsurf_save       _NULL
 end type cosp_from_rad_block_type

!--- Cosp Rad control type
 public cosp_from_rad_control_type
 type cosp_from_rad_control_type
     logical :: step_to_call_cosp
 end type cosp_from_rad_control_type

!--- Cosp_rad type definition
 public cosp_from_rad_type
 type cosp_from_rad_type
     type (cosp_from_rad_control_type) :: control
     type (cosp_from_rad_block_type), dimension(:), _ALLOCATABLE :: block _NULL
 end type cosp_from_rad_type

public :: exch_rad_phys_state, &
          alloc_clouds_from_moist_type, dealloc_clouds_from_moist_type, &
          alloc_cosp_from_rad_type, dealloc_cosp_from_rad_type, &
          alloc_radiation_flux_type, dealloc_radiation_flux_type

contains


 subroutine alloc_cosp_from_rad_type (Cosp_rad, Exch_ctrl, Atm_block)
   type (cosp_from_rad_type),    intent(inout) :: Cosp_rad(:)
   type (exchange_control_type), intent(in)    :: Exch_ctrl
   type (block_control_type),    intent(in)    :: Atm_block
!--- local variables
   integer :: n, nb, ix, jx, npz

     npz = Atm_block%npz
     ncol = Exch_ctrl%ncol
     do n = 1, size(Cosp_rad,1)
       allocate (Cosp_rad(n)%block(Atm_block%nblks))
       if (Exch_ctrl%do_cosp .or. Exch_ctrl%do_modis_yim) then
         do nb = 1,Atm_block%nblks
           ix = Atm_block%ibe(nb)-Atm_block%ibs(nb)+1
           jx = Atm_block%jbe(nb)-Atm_block%jbs(nb)+1
           allocate ( Cosp_rad(n)%block(nb)%tau_stoch       (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%lwem_stoch      (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_cloud_type(ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_conc_drop (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_conc_ice  (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_size_drop (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_size_ice  (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%mr_ozone        (ix,jx,npz),      &
                      Cosp_rad(n)%block(nb)%tsurf_save      (ix,jx) ,         &
                      Cosp_rad(n)%block(nb)%daytime         (ix,jx)          )
      ! initial values
           Cosp_rad(n)%block(nb)%tau_stoch        = 0.
           Cosp_rad(n)%block(nb)%lwem_stoch       = 0.
           Cosp_rad(n)%block(nb)%stoch_cloud_type = 0.
           Cosp_rad(n)%block(nb)%stoch_conc_drop  = 0.
           Cosp_rad(n)%block(nb)%stoch_conc_ice   = 0.
           Cosp_rad(n)%block(nb)%stoch_size_drop  = 0.
           Cosp_rad(n)%block(nb)%stoch_size_ice   = 0.
           Cosp_rad(n)%block(nb)%mr_ozone         = 0.
           Cosp_rad(n)%block(nb)%tsurf_save       = 0.
           Cosp_rad(n)%block(nb)%daytime          = 0.
         end do
       endif
     end do

 end subroutine alloc_cosp_from_rad_type



 subroutine dealloc_cosp_from_rad_type (Cosp_rad, Exch_ctrl)
   type (cosp_from_rad_type), intent(inout) :: Cosp_rad(:)
   type (exchange_control_type), intent(in) :: Exch_ctrl
!--- local variables
   integer :: n, nb

    do n = 1, size(Cosp_rad,1)
      if (Exch_ctrl%do_cosp .or. Exch_ctrl%do_modis_yim) then
        do nb = 1, size(Cosp_rad(n)%block,1)
          deallocate (Cosp_rad(n)%block(nb)%tau_stoch,        &
                      Cosp_rad(n)%block(nb)%lwem_stoch,       &
                      Cosp_rad(n)%block(nb)%stoch_cloud_type, &
                      Cosp_rad(n)%block(nb)%stoch_conc_drop,  &
                      Cosp_rad(n)%block(nb)%stoch_conc_ice,   &
                      Cosp_rad(n)%block(nb)%stoch_size_drop,  &
                      Cosp_rad(n)%block(nb)%stoch_size_ice,   &
                      Cosp_rad(n)%block(nb)%mr_ozone,         &
                      Cosp_rad(n)%block(nb)%tsurf_save,       &
                      Cosp_rad(n)%block(nb)%daytime)
        end do
      endif
      deallocate(Cosp_rad(n)%block)
    end do

 end subroutine dealloc_cosp_from_rad_type



 subroutine alloc_radiation_flux_type (Rad_flux, nonzero_init, Atm_block)
   type (radiation_flux_type), intent(inout) :: Rad_flux(:)
   logical,                    intent(in)    :: nonzero_init
   type (block_control_type),  intent(in)    :: Atm_block
!--- local variables
   integer :: n, nb
   integer :: ix, jx, npz

    npz = Atm_block%npz
    do n = 1, size(Rad_flux,1)
      allocate (Rad_flux(n)%block(Atm_block%nblks))
      do nb = 1, Atm_block%nblks
        ix = Atm_block%ibe(nb) - Atm_block%ibs(nb) + 1
        jx = Atm_block%jbe(nb) - Atm_block%jbs(nb) + 1
        call Rad_flux(n)%block(nb)%alloc ( ix, jx, npz, nonzero_init )
      end do
    end do

 end subroutine alloc_radiation_flux_type

!######################################################################

 subroutine dealloc_radiation_flux_type (Rad_flux)
   type (radiation_flux_type), intent(inout) :: Rad_flux(:)
!--- local variables
   integer :: n, nb
!---------------------------------------------------------------------
!    deallocate the variables
!---------------------------------------------------------------------
      do n = 1, size(Rad_flux,1)
        do nb = 1, size(Rad_flux(n)%block,1)
          call Rad_flux(n)%block(nb)%dealloc
        enddo
        deallocate (Rad_flux(n)%block)
      end do
 end subroutine dealloc_radiation_flux_type

!######################################################################

 subroutine alloc_radiation_flux_block_type (Rad_flux_block, ix, jx, kx, nonzero_init)
   class(radiation_flux_block_type), intent(inout) :: Rad_flux_block
   integer,                          intent(in)    :: ix, jx, kx
   logical, optional,                intent(in)    :: nonzero_init
   logical :: nonzero_init_local

   nonzero_init_local = .false.
   if (present(nonzero_init)) nonzero_init_local = nonzero_init

      allocate (Rad_flux_block%tdt_rad               (ix,jx,kx), &
                Rad_flux_block%tdt_lw                (ix,jx,kx), &
                Rad_flux_block%flux_sw               (ix,jx),    &
                Rad_flux_block%flux_sw_dir           (ix,jx),    &
                Rad_flux_block%flux_sw_dif           (ix,jx),    &
                Rad_flux_block%flux_sw_down_vis_dir  (ix,jx),    &
                Rad_flux_block%flux_sw_down_vis_dif  (ix,jx),    &
                Rad_flux_block%flux_sw_down_total_dir(ix,jx),    &
                Rad_flux_block%flux_sw_down_total_dif(ix,jx),    &
                Rad_flux_block%flux_sw_vis           (ix,jx),    &
                Rad_flux_block%flux_sw_vis_dir       (ix,jx),    &
                Rad_flux_block%flux_sw_vis_dif       (ix,jx),    &
                Rad_flux_block%flux_lw               (ix,jx),    &
                Rad_flux_block%coszen                (ix,jx),    &
                Rad_flux_block%extinction            (ix,jx,kx)  )

      ! initial values - only used at t=0 when there is no restart
      if (nonzero_init_local) then
          Rad_flux_block%tdt_rad                = 0.0
          Rad_flux_block%tdt_lw                 = 0.0
          Rad_flux_block%flux_sw                = 50.0
          Rad_flux_block%flux_sw_dir            = 25.0
          Rad_flux_block%flux_sw_dif            = 25.0
          Rad_flux_block%flux_sw_down_vis_dir   = 12.5 * 1.1
          Rad_flux_block%flux_sw_down_vis_dif   = 12.5 * 1.1
          Rad_flux_block%flux_sw_down_total_dir = 25.0 * 1.1
          Rad_flux_block%flux_sw_down_total_dif = 25.0 * 1.1
          Rad_flux_block%flux_sw_vis            = 25.0
          Rad_flux_block%flux_sw_vis_dir        = 12.5
          Rad_flux_block%flux_sw_vis_dif        = 12.5
          Rad_flux_block%flux_lw                = 50.0
          Rad_flux_block%coszen                 = 0.50
          Rad_flux_block%extinction             = 0.0
      else
          Rad_flux_block%tdt_rad                = 0.0
          Rad_flux_block%tdt_lw                 = 0.0
          Rad_flux_block%flux_sw                = 0.0 !50.0
          Rad_flux_block%flux_sw_dir            = 0.0 !25.0
          Rad_flux_block%flux_sw_dif            = 0.0 !25.0
          Rad_flux_block%flux_sw_down_vis_dir   = 0.0 !12.5 * 1.1
          Rad_flux_block%flux_sw_down_vis_dif   = 0.0 !12.5 * 1.1
          Rad_flux_block%flux_sw_down_total_dir = 0.0 !25.0 * 1.1
          Rad_flux_block%flux_sw_down_total_dif = 0.0 !25.0 * 1.1
          Rad_flux_block%flux_sw_vis            = 0.0 !25.0
          Rad_flux_block%flux_sw_vis_dir        = 0.0 !12.5
          Rad_flux_block%flux_sw_vis_dif        = 0.0 !12.5
          Rad_flux_block%flux_lw                = 0.0 !50.0
          Rad_flux_block%coszen                 = 0.0 !0.50
          Rad_flux_block%extinction             = 0.0 !0.0
      endif

 end subroutine alloc_radiation_flux_block_type

!######################################################################

 subroutine dealloc_radiation_flux_block_type (Rad_flux_block)
   class(radiation_flux_block_type), intent(inout) :: Rad_flux_block

      deallocate (Rad_flux_block%tdt_rad,                &
                  Rad_flux_block%tdt_lw,                 &
                  Rad_flux_block%flux_sw,                &
                  Rad_flux_block%flux_sw_dir,            &
                  Rad_flux_block%flux_sw_dif,            &
                  Rad_flux_block%flux_sw_down_vis_dir,   &
                  Rad_flux_block%flux_sw_down_vis_dif,   &
                  Rad_flux_block%flux_sw_down_total_dir, &
                  Rad_flux_block%flux_sw_down_total_dif, &
                  Rad_flux_block%flux_sw_vis,            &
                  Rad_flux_block%flux_sw_vis_dir,        &
                  Rad_flux_block%flux_sw_vis_dif,        &
                  Rad_flux_block%flux_lw,                &
                  Rad_flux_block%coszen,                 &
                  Rad_flux_block%extinction              )

 end subroutine dealloc_radiation_flux_block_type

!######################################################################

 subroutine alloc_clouds_from_moist_type (Moist_clouds, Exch_ctrl, Atm_block)
   type (clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
   type (exchange_control_type),  intent(in)    :: Exch_ctrl
   type (block_control_type),     intent(in)    :: Atm_block
!--- local variables
   integer :: n, nb, npz
   integer :: ix, jx
   integer :: nc

   npz = Atm_block%npz
!-----------------------------------------------------------------------
!    allocate derived-type that stores cloud properties
!    return from moist processes
!-----------------------------------------------------------------------
    do n=1,size(Moist_clouds,1)
      allocate(Moist_clouds(n)%block(Atm_block%nblks))
      do nb = 1, Atm_block%nblks
        ix = Atm_block%ibe(nb) - Atm_block%ibs(nb) + 1
        jx = Atm_block%jbe(nb) - Atm_block%jbs(nb) + 1

        allocate(Moist_clouds(n)%block(nb)%Cloud_data(Exch_ctrl%ncld)) ! may want to compute locally from flag in Exch_ctrl
        nc = 0

!-------------------------------------------------------------------------
!    Initialize here so that these may be variables may be referenced even 
!    when scheme not active.
!-------------------------------------------------------------------------
        Moist_clouds(n)%block(nb)%index_strat = 0
        Moist_clouds(n)%block(nb)%index_donner_cell = 0
        Moist_clouds(n)%block(nb)%index_donner_meso = 0
        Moist_clouds(n)%block(nb)%index_uw_conv     = 0

        if (Exch_ctrl%doing_prog_clouds) then
           nc = nc+1
           Moist_clouds(n)%block(nb)%index_strat = nc
           call alloc_cloud_scheme_data_type('strat_cloud',ix,jx,npz,Moist_clouds(n)%block(nb)%Cloud_data(nc))
        endif

        if (Exch_ctrl%doing_donner) then
           ! cell
           nc = nc+1
           Moist_clouds(n)%block(nb)%index_donner_cell = nc
           call alloc_cloud_scheme_data_type('donner_cell',ix,jx,npz,Moist_clouds(n)%block(nb)%Cloud_data(nc))
           ! meso
           nc = nc+1
           Moist_clouds(n)%block(nb)%index_donner_meso = nc
           call alloc_cloud_scheme_data_type('donner_meso',ix,jx,npz,Moist_clouds(n)%block(nb)%Cloud_data(nc))
        endif

        if (Exch_ctrl%doing_uw_conv) then
           nc = nc+1
           Moist_clouds(n)%block(nb)%index_uw_conv = nc
           call alloc_cloud_scheme_data_type('uw_conv',ix,jx,npz,Moist_clouds(n)%block(nb)%Cloud_data(nc))
        endif
      enddo
    enddo

 end subroutine alloc_clouds_from_moist_type

!######################################################################

 subroutine alloc_cloud_scheme_data_type (scheme, id, jd, kd, Cloud_data )
   character(len=*),              intent(in)    :: scheme
   integer,                       intent(in)    :: id, jd, kd
   type (cloud_scheme_data_type), intent(inout) :: Cloud_data

   logical :: done_allocation = .false.

!-----------------------------------------------------------------------
!    allocate derived-type that stores cloud properties
!    return from moist processes
!-----------------------------------------------------------------------

     ! properties common to all cloud schemes
      allocate (Cloud_data%cloud_area     (id, jd, kd) )
      allocate (Cloud_data%liquid_amt     (id, jd, kd) )
      allocate (Cloud_data%ice_amt        (id, jd, kd) )
      allocate (Cloud_data%droplet_number (id, jd, kd) )
     !Cloud_data%cloud_area     = 0.0
     !Cloud_data%liquid_amt     = 0.0
     !Cloud_data%ice_amt        = 0.0
     !Cloud_data%droplet_number = 0.0

     ! properties specific to large-scale/stratiform clouds
      if (lowercase(trim(scheme)) .eq. 'strat_cloud') then

          allocate (Cloud_data%ice_number     (id, jd, kd) )
          allocate (Cloud_data%rain           (id, jd, kd) )
          allocate (Cloud_data%snow           (id, jd, kd) )
          allocate (Cloud_data%rain_size      (id, jd, kd) )
          allocate (Cloud_data%snow_size      (id, jd, kd) )
          Cloud_data%cloud_area      = -99.
          Cloud_data%liquid_amt      = -99.
          Cloud_data%ice_amt         = -99.
          Cloud_data%droplet_number  = -99.
          Cloud_data%ice_number = -99.
          Cloud_data%rain = 0.
          Cloud_data%snow = 0.
          Cloud_data%rain_size = 0.
          Cloud_data%snow_size = 0.
          Cloud_data%scheme_name = lowercase(trim(scheme))
          done_allocation = .true.
      else
          Cloud_data%cloud_area      = 0.
          Cloud_data%liquid_amt      = 0.
          Cloud_data%ice_amt         = 0.
          Cloud_data%droplet_number  = 0.
      end if

      ! properties specific to donner deep clouds (both cell and meso)
      if (lowercase(trim(scheme)) .eq. 'donner_cell' .or. &
          lowercase(trim(scheme)) .eq. 'donner_meso') then

          allocate (Cloud_data%liquid_size (id, jd, kd) )
          allocate (Cloud_data%ice_size    (id, jd, kd) )
          allocate (Cloud_data%nsum_out    (id, jd) )
          Cloud_data%liquid_size = 0.
          Cloud_data%ice_size    = 0.
          Cloud_data%nsum_out    = 1
          Cloud_data%scheme_name = lowercase(trim(scheme))
          done_allocation = .true.
      endif

      ! properties specific to uw shallow convective clouds
      if (lowercase(trim(scheme)) .eq. 'uw_conv') then

          allocate (Cloud_data%ice_number (id, jd, kd) )
          Cloud_data%ice_number = 0.
          Cloud_data%scheme_name = lowercase(trim(scheme))
          done_allocation = .true.
      endif

      ! verify that the arrays were allocated
      if (.not.done_allocation) then
          call error_mesg ('physics_radiation_exch_mod', &
                           'invalid cloud scheme name', FATAL)
      endif

!----------------------------------------------------------------------

end subroutine alloc_cloud_scheme_data_type

!######################################################################

 subroutine dealloc_clouds_from_moist_type (Moist_clouds, Exch_ctrl)
   type (clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
   type (exchange_control_type),  intent(in)    :: Exch_ctrl
!--- local variables
   integer :: n, nb, nc
!--------------------------------------------------------------------
!    deallocate variables
!--------------------------------------------------------------------
    do n=1,size(Moist_clouds,1)
      do nb = 1, size(Moist_clouds(n)%block,1)
        do nc = 1, size(Moist_clouds(n)%block(nb)%Cloud_data,1)

          ! deallocate arrays common to all cloud schemes
          deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%cloud_area    )
          deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%liquid_amt    )
          deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%ice_amt       )
          deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%droplet_number)

          ! properties specific to large-scale/stratiform clouds
          if (trim(Moist_clouds(n)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'strat_cloud') then
            deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%ice_number)
            deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%rain      )
            deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%snow      )
            deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%rain_size )
            deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%snow_size )
          endif

          ! properties specific to donner deep clouds (both cell and meso)
          if (trim(Moist_clouds(n)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'donner_cell' .or. &
              trim(Moist_clouds(n)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'donner_meso') then
            deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%liquid_size)
            deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%ice_size   )
            deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%nsum_out   )
          endif

          ! properties specific to uw shallow convective clouds
          if (trim(Moist_clouds(n)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'uw_conv') then
            deallocate (Moist_clouds(n)%block(nb)%Cloud_data(nc)%ice_number)
          endif

          Moist_clouds(n)%block(nb)%Cloud_data(nc)%scheme_name = ' '

        enddo
        deallocate (Moist_clouds(n)%block(nb)%Cloud_data)
      enddo
      deallocate (Moist_clouds(n)%block)
    enddo

 end subroutine dealloc_clouds_from_moist_type

!######################################################################

! update radiation flux states
 subroutine exch_rad_phys_state (Moist_clouds, Cosp_rad, Rad_flux, Exch_ctrl)
   type (clouds_from_moist_type), intent(inout) :: Moist_clouds(2)
   type (cosp_from_rad_type),     intent(inout) :: Cosp_rad(2)
   type (radiation_flux_type),    intent(inout) :: Rad_flux(2)
   type (exchange_control_type),  intent(in)    :: Exch_ctrl
!--- local variables
   integer :: nb

    Cosp_rad(1)%control%step_to_call_cosp = Cosp_rad(2)%control%step_to_call_cosp
    Rad_flux(1)%control%do_rad = Rad_flux(2)%control%do_rad
!$OMP parallel do default(shared) private(nb)
    do nb = 1, size(Rad_flux(1)%block)
       Rad_flux(1)%block(nb)%tdt_rad                = Rad_flux(2)%block(nb)%tdt_rad
       Rad_flux(1)%block(nb)%tdt_lw                 = Rad_flux(2)%block(nb)%tdt_lw
       Rad_flux(1)%block(nb)%flux_sw                = Rad_flux(2)%block(nb)%flux_sw
       Rad_flux(1)%block(nb)%flux_sw_dir            = Rad_flux(2)%block(nb)%flux_sw_dir
       Rad_flux(1)%block(nb)%flux_sw_dif            = Rad_flux(2)%block(nb)%flux_sw_dif
       Rad_flux(1)%block(nb)%flux_sw_down_vis_dir   = Rad_flux(2)%block(nb)%flux_sw_down_vis_dir
       Rad_flux(1)%block(nb)%flux_sw_down_vis_dif   = Rad_flux(2)%block(nb)%flux_sw_down_vis_dif
       Rad_flux(1)%block(nb)%flux_sw_down_total_dir = Rad_flux(2)%block(nb)%flux_sw_down_total_dir
       Rad_flux(1)%block(nb)%flux_sw_down_total_dif = Rad_flux(2)%block(nb)%flux_sw_down_total_dif
       Rad_flux(1)%block(nb)%flux_sw_vis            = Rad_flux(2)%block(nb)%flux_sw_vis
       Rad_flux(1)%block(nb)%flux_sw_vis_dir        = Rad_flux(2)%block(nb)%flux_sw_vis_dir
       Rad_flux(1)%block(nb)%flux_sw_vis_dif        = Rad_flux(2)%block(nb)%flux_sw_vis_dif
       Rad_flux(1)%block(nb)%flux_lw                = Rad_flux(2)%block(nb)%flux_lw
       Rad_flux(1)%block(nb)%coszen                 = Rad_flux(2)%block(nb)%coszen
       Rad_flux(1)%block(nb)%extinction             = Rad_flux(2)%block(nb)%extinction
       if (Exch_ctrl%do_cosp .or. Exch_ctrl%do_modis_yim) then
         Cosp_rad(1)%block(nb)%tau_stoch        = Cosp_rad(2)%block(nb)%tau_stoch
         Cosp_rad(1)%block(nb)%lwem_stoch       = Cosp_rad(2)%block(nb)%lwem_stoch
         Cosp_rad(1)%block(nb)%stoch_cloud_type = Cosp_rad(2)%block(nb)%stoch_cloud_type
         Cosp_rad(1)%block(nb)%stoch_conc_drop  = Cosp_rad(2)%block(nb)%stoch_conc_drop
         Cosp_rad(1)%block(nb)%stoch_conc_ice   = Cosp_rad(2)%block(nb)%stoch_conc_ice
         Cosp_rad(1)%block(nb)%stoch_size_drop  = Cosp_rad(2)%block(nb)%stoch_size_drop
         Cosp_rad(1)%block(nb)%stoch_size_ice   = Cosp_rad(2)%block(nb)%stoch_size_ice 
         Cosp_rad(1)%block(nb)%mr_ozone         = Cosp_rad(2)%block(nb)%mr_ozone
         Cosp_rad(1)%block(nb)%daytime          = Cosp_rad(2)%block(nb)%daytime 
         Cosp_rad(1)%block(nb)%tsurf_save       = Cosp_rad(2)%block(nb)%tsurf_save
       endif

       do nc = 1, size(Moist_clouds(1)%block(nb)%Cloud_data,1)

         ! common to all cloud schemes
         Moist_clouds(2)%block(nb)%Cloud_data(nc)%cloud_area     = Moist_clouds(1)%block(nb)%Cloud_data(nc)%cloud_area
         Moist_clouds(2)%block(nb)%Cloud_data(nc)%liquid_amt     = Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_amt
         Moist_clouds(2)%block(nb)%Cloud_data(nc)%ice_amt        = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_amt
         Moist_clouds(2)%block(nb)%Cloud_data(nc)%droplet_number = Moist_clouds(1)%block(nb)%Cloud_data(nc)%droplet_number
         Moist_clouds(2)%block(nb)%Cloud_data(nc)%scheme_name    = Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name

         ! properties specific to large-scale/stratiform clouds
         if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'strat_cloud') then
           Moist_clouds(2)%block(nb)%Cloud_data(nc)%ice_number = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number
           Moist_clouds(2)%block(nb)%Cloud_data(nc)%rain       = Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain
           Moist_clouds(2)%block(nb)%Cloud_data(nc)%snow       = Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow
           Moist_clouds(2)%block(nb)%Cloud_data(nc)%rain_size  = Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain_size
           Moist_clouds(2)%block(nb)%Cloud_data(nc)%snow_size  = Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow_size
         endif

         ! properties specific to donner deep clouds (both cell and meso)
         if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'donner_cell' .or. &
             trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'donner_meso') then
           Moist_clouds(2)%block(nb)%Cloud_data(nc)%liquid_size = Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_size
           Moist_clouds(2)%block(nb)%Cloud_data(nc)%ice_size    = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_size
           Moist_clouds(2)%block(nb)%Cloud_data(nc)%nsum_out    = Moist_clouds(1)%block(nb)%Cloud_data(nc)%nsum_out
         endif

         ! properties specific to uw shallow convective clouds
         if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'uw_conv') then
           Moist_clouds(2)%block(nb)%Cloud_data(nc)%ice_number = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number
         endif

       enddo
     enddo

 end subroutine exch_rad_phys_state

end module physics_radiation_exch_mod
