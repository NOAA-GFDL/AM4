
module cloudrad_types_mod

!--------------------------------------------------------------------

implicit none 
private 

!--------------------------------------------------------------------
!---- public data structures ----
!--------------------------------------------------------------------

public cld_specification_type

!
!  The following elements are components of the cld_specification_type:
!
!    %cmxolw          amount of maximally overlapped longwave clouds [ dimensionless ]
!    %crndlw          amount of randomly overlapped longwave clouds [ dimensionless ]
!    %nmxolw          number of maximally overlapped longwave clouds
!                     in each grid column
!    %nrndlw          number of maximally overlapped longwave clouds
!                     in each grid column
!    %camtsw          shortwave cloud amount. the sum of the maximally
!                     overlapped and randomly overlapped longwave cloud
!                     amounts [ dimensionless ]
!    %ncldsw          number of shortwave clouds in each grid column
!    %camtsw_band     shortwave cloud amount. the sum of the maximally
!                     overlapped and randomly overlapped longwave cloud amounts,
!                     differing with sw parameterization band [ dimensionless ]
!    %crndlw_band     amount of randomly overlapped longwave clouds,
!                     differing with lw parameterization band [ dimensionless ]
!    %hi_cloud        logical mask for high clouds 
!    %mid_cloud       logical mask for middle clouds
!    %low_cloud       logical mask for low clouds
!    %ice_cloud       logical mask for ice clouds
!    %iwp             ice water path  [ kg / m**2 ]
!    %lwp             liquid water path [ kg / m**2 ]
!    %reff_liq        effective cloud drop radius  used with
!                     bulk cloud physics scheme [ microns ]
!    %reff_ice        effective ice crystal radius used with
!                     bulk cloud physics scheme [ microns ]
!    %reff_liq_micro  effective cloud drop radius used with 
!                     microphysically based scheme [ microns ]
!    %reff_ice_micro  effective ice crystal radius used with
!                     microphysically based scheme [ microns ]
!    %tau             extinction optical path  [ dimensionless ]
!    %liq_frac        fraction of cloud in a box which is liquid [ dimensionless ]
!    %cld_thickness   number of model layers contained in cloud  
!    %cloud_water     liquid cloud content [ kg liq / kg air ]
!    %cloud_ice       ice cloud content [ kg ice / kg air ]
!    %cloud_area      saturated volume fraction [ dimensionless ]

type cld_specification_type
!BWreal, dimension(:,:,:,:),  pointer :: tau=>NULL(),  &
   real, dimension(:,:,:,:),  pointer :: camtsw_band=>NULL(), &
                                         crndlw_band=>NULL(), &
                                         lwp_lw_band=>NULL(), &
                                         iwp_lw_band=>NULL(), &
                                         lwp_sw_band=>NULL(), &
                                         iwp_sw_band=>NULL(), &
                                         reff_liq_lw_band=>NULL(),   &
                                         reff_ice_lw_band=>NULL(), &
                                         reff_liq_sw_band=>NULL(),   &
                                         reff_ice_sw_band=>NULL()
   real, dimension(:,:,:),    pointer :: lwp=>NULL(),   &
                                         iwp=>NULL(),  &
                                         reff_liq=>NULL(),   &
                                         reff_ice=>NULL(), &
                                         reff_liq_lim=>NULL(),   &
                                         reff_ice_lim=>NULL(), &
                                         liq_frac=>NULL(), &
                                         cloud_water=>NULL(), &
                                         cloud_ice=>NULL(),  &
                                         cloud_area=>NULL(), &
                                         cloud_droplet=>NULL(), &
                                         cloud_ice_num=>NULL(), &
                                         rain =>NULL(), &
                                         snow =>NULL(), &
                                         rain_size =>NULL(), &
                                         snow_size =>NULL(), &
                                         reff_liq_micro=>NULL(),   &
                                         reff_ice_micro=>NULL(),&
                                         camtsw=>NULL(),   &
                                         cmxolw=>NULL(),  &
                                         crndlw=>NULL()
   integer, dimension(:,:,:), pointer :: cld_thickness=>NULL()
   integer, dimension(:,:,:,:), pointer :: stoch_cloud_type=>NULL()
   integer, dimension(:,:,:,:), pointer :: cld_thickness_lw_band=>NULL()
   integer, dimension(:,:,:,:), pointer :: cld_thickness_sw_band=>NULL()
   integer, dimension(:,:),   pointer :: ncldsw=>NULL(),   &
                                         nmxolw=>NULL(),&
                                         nrndlw=>NULL()
   integer, dimension(:,:,:), pointer :: ncldsw_band=>NULL(),   &
                                         nrndlw_band=>NULL()
   logical, dimension(:,:,:), pointer :: hi_cloud=>NULL(),   &
                                         mid_cloud=>NULL(),  &
                                         low_cloud=>NULL(),   &
                                         ice_cloud=>NULL()

   contains
      procedure :: alloc=>cloud_spec_alloc
      procedure :: dealloc=>cloud_spec_dealloc
end type cld_specification_type

!--------------------------------------------------------------------

public cldrad_properties_type

!
!  The components of the cldrad_properties_type structure are:
!
!     %emmxolw    longwave cloud emissivity for maximally overlapped clouds
!                 [ dimensionless ] 
!     %emrndlw    longwave cloud emissivity for randomly overlapped clouds
!                 [ dimensionless ]
!     %abscoeff   combined absorption coefficient for clouds in each of the
!                 longwave frequency bands [ km**(-1) ]
!
!     %cldext     parameterization band values of the cloud extinction coefficient
!                 [ km**(-1) ]   
!     %cldsct     parameterization band values of the cloud scattering coefficient
!                 [ km**(-1) ]
!     %cldasymm   parameterization band values of the asymmetry factor
!                 [ dimensionless ]
!
!  The following are no longer used:
!     %cldemiss   longwave emissivity calculated using abscoeff
!                 [ dimensionless ]
!     %cirabsw    absorptivity of clouds in the infrared frequency band
!                 may be zenith angle dependent
!                 [ dimensionless ]
!     %cirrfsw    reflectivity of clouds in the infrared frequency band
!                 may be zenith angle dependent
!                 [ dimensionless ]
!     %cvisrfsw   reflectivity of clouds in the visible frequency band
!                 may be zenith angle dependent
!                 [ dimensionless ]

type cldrad_properties_type
     real, dimension(:,:,:,:,:), pointer :: cldext=>NULL(),   &
                                            cldasymm=>NULL(), &
                                            cldsct=>NULL()
     real, dimension(:,:,:,:,:), pointer :: emmxolw=>NULL(),  &
                                            emrndlw=>NULL(),  &
                                            abscoeff=>NULL(), &
                                            cldemiss=>NULL()
     real, dimension(:,:,:),     pointer :: cirabsw=>NULL(), &
                                            cirrfsw=>NULL(), &
                                            cvisrfsw=>NULL()

     contains
       procedure :: alloc=>cldrad_properties_alloc
       procedure :: dealloc=>cldrad_properties_dealloc
end type cldrad_properties_type

!--------------------------------------------------------------------

public cloudrad_control_type

type cloudrad_control_type
    logical :: do_pred_cld_microphys
    logical :: do_presc_cld_microphys
    logical :: do_bulk_microphys
    logical :: do_sw_micro
    logical :: do_lw_micro
    logical :: do_strat_clouds
    logical :: do_no_clouds
    logical :: do_donner_deep_clouds
    logical :: do_uw_clouds
    logical :: do_random_overlap
    logical :: do_max_random_overlap
    logical :: do_stochastic_clouds
    logical :: use_temp_for_seed
    logical :: do_ica_calcs
    logical :: do_liq_num
    logical :: do_ice_num
    logical :: using_fu2007
    integer :: num_sw_cloud_bands
    integer :: num_lw_cloud_bands
end type cloudrad_control_type

!--------------------------------------------------------------------

public microphysics_type

!
!   The following elements are components of the microphysics_type:
!
!     %conc_ice     ice particle concentration [ g / m**3 ]
!     %conc_drop    cloud droplet concentration [ g / m**3 ]
!     %conc_rain    rain drop concentration [ g / m**3 ]
!     %conc_snow    snow concentration [ g / m**3 ]
!     %size_ice     effective ice crystal diameter [ microns ]
!     %size_drop    effective cloud drop diameter [ microns ]
!     %size_rain    effective rain drop diameter [ microns ]
!     %size_snow    effective snow flake diameter [ microns ]
!     %cldamt       total cloud fraction (crystal + droplet)
!                   [ dimensionless ]
!     %lw_stoch_conc_ice
!                   ice particle concentration as a function of
!                   lw parameterization band [ g / m**3 ]
!     %lw_stoch_conc_drop
!                   cloud droplet concentration as a function of
!                   lw parameterization band [ g / m**3 ]
!     %lw_stoch_size_ice
!                   effective ice crystal diameter as a function
!                   of lw parameterization band [ microns ]
!     %lw_stoch_size_drop
!                   effective cloud drop diameter as a function of
!                   lw parameterization band [ microns ]
!     %lw_stoch_cldamt
!                   total cloud fraction (crystal + droplet) as a
!                   function of lw parameterization band [ dimensionless ]
!     %sw_stoch_conc_ice
!                   ice particle concentration as a function of
!                   sw parameterization band [ g / m**3 ]
!     %sw_stoch_conc_drop
!                   cloud droplet concentration as a function of
!                   sw parameterization band [ g / m**3 ]
!     %sw_stoch_size_ice
!                   effective ice crystal diameter as a function
!                   of sw parameterization band [ microns ]
!     %sw_stoch_size_drop
!                   effective cloud drop diameter as a function of
!                   sw parameterization band [ microns ]
!     %sw_stoch_cldamt
!                   total cloud fraction (crystal + droplet) as a
!                   function of sw parameterization band [ dimensionless ]

type microphysics_type
character(len=64) :: scheme_name
logical :: use_for_diag
real, dimension(:,:,:), pointer    :: conc_ice=>NULL(),   &
                                      conc_drop=>NULL(),      &
                                      size_ice=>NULL(),   &
                                      size_drop=>NULL(),     &
                                      size_snow=>NULL(),   &
                                      conc_snow=>NULL(),     &
                                      size_rain=>NULL(),     &
                                      conc_rain=>NULL(),   &
                                      cldamt=>NULL(),      &
                                      droplet_number=>NULL(), &
                                      ice_number=>NULL()
!  The following are activated for stochastic clouds
real, dimension(:,:,:,:), pointer :: stoch_conc_ice=>NULL(),   &
                                     stoch_conc_drop=>NULL(),  &
                                     stoch_size_ice=>NULL(),   &
                                     stoch_size_drop=>NULL(),  &
                                     stoch_cldamt=>NULL(),     &
                                     stoch_droplet_number=>NULL(), &
                                     stoch_ice_number=>NULL()
integer, dimension(:,:,:,:), pointer ::  stoch_cloud_type=>NULL()

!  In practice, we allocate a single set of columns for the
!  stochastic clouds, then point to sections of the larger array
!  with the lw_ and sw_ pointer arrays. 
!  i.e., lw_stoch_conc_ice => stoch_conc_ice(:, :, :, 1:numLwBands)

real, dimension(:,:,:,:), pointer :: lw_stoch_conc_ice=>NULL(),   &
                                     lw_stoch_conc_drop=>NULL(),  &
                                     lw_stoch_size_ice=>NULL(),   &
                                     lw_stoch_size_drop=>NULL(),  &
                                     lw_stoch_cldamt=>NULL(),     &
                                     lw_stoch_droplet_number=>NULL(), &
                                     lw_stoch_ice_number=>NULL(), &
                                     sw_stoch_conc_ice=>NULL(),   &
                                     sw_stoch_conc_drop=>NULL(),  &
                                     sw_stoch_size_ice=>NULL(),   &
                                     sw_stoch_size_drop=>NULL(),  &
                                     sw_stoch_cldamt=>NULL(),     &
                                     sw_stoch_droplet_number=>NULL(), &
                                     sw_stoch_ice_number=>NULL()

!  The following are only activated for large-scale cloud diagnostics

real, dimension(:,:,:), pointer    :: lsc_conc_drop=>NULL(),      &
                                      lsc_size_drop=>NULL(),     &
                                      lsc_cldamt=>NULL(),      &
                                      lsc_droplet_number=>NULL()

   contains
      procedure :: microphysics_alloc
      procedure :: microphysics_alloc_diag
      generic :: alloc=>microphysics_alloc, microphysics_alloc_diag
      procedure :: dealloc=>microphysics_dealloc
end type microphysics_type

!--------------------------------------------------------------------

public microrad_properties_type

!
!  The components of the microrad_structure are:
!
!     %cldext    parameterization band values of the cloud extinction coefficient
!                [ km**(-1) ]   
!     %cldsct    parameterization band values of the cloud scattering coefficient
!                [ km**(-1) ]
!     %cldasymm  parameterization band values of the asymmetry factor
!                [ dimensionless ]
!
!     %abscoeff  combined absorption coefficient for clouds in each of the
!                longwave frequency bands [ km**(-1) ]

type microrad_properties_type
   character(len=64) :: scheme_name
   real, dimension(:,:,:,:), pointer :: cldext=>NULL(),  &
                                        cldsct=>NULL(), &
                                        cldasymm=>NULL(),    &
                                        abscoeff=>NULL()

   contains
    ! procedure :: alloc=>microrad_properties_alloc
      procedure :: dealloc=>microrad_properties_dealloc
end type microrad_properties_type

public :: microrad_properties_alloc
!--------------------------------------------------------------------

CONTAINS

!####################################################################
!####################################################################
!             constructor/destructor routines
!####################################################################
!####################################################################

subroutine cloud_spec_alloc (Cld_spec, ix, jx, kx, Cldrad_control)

class(cld_specification_type), intent(inout) :: Cld_spec
integer,                       intent(in)    :: ix, jx, kx
type(cloudrad_control_type),   intent(in)    :: Cldrad_control

!--------------------------------------------------------------------
! local variables

integer :: n, nswcb, nlwcb

!---------------------------------------------------------------------

      nswcb = Cldrad_control%num_sw_cloud_bands
      nlwcb = Cldrad_control%num_lw_cloud_bands

!---------------------------------------------------------------------
!    allocate arrays to hold the cloud fractions seen by the shortwave
!    and the random and maximum overlap fractions seen by the longwave
!    radiation, and then the number of each of these types of cloud in
!    each column. initialize the cloud fractions and number of clouds
!    to zero.
!---------------------------------------------------------------------
      allocate ( Cld_spec%camtsw (ix, jx, kx ) )
      allocate ( Cld_spec%cmxolw (ix, jx, kx ) )
      allocate ( Cld_spec%crndlw (ix, jx, kx ) )
      allocate ( Cld_spec%ncldsw (ix, jx     ) )  
      allocate ( Cld_spec%nmxolw (ix, jx     ) )  
      allocate ( Cld_spec%nrndlw (ix, jx     ) )  
      Cld_spec%cmxolw(:,:,:) = 0.0
      Cld_spec%crndlw(:,:,:) = 0.0
      Cld_spec%camtsw(:,:,:) = 0.0
      Cld_spec%nmxolw (:,:)  = 0  
      Cld_spec%nrndlw (:,:)  = 0  
      Cld_spec%ncldsw (:,:)  = 0  

      if (Cldrad_control%do_stochastic_clouds) then 
        ! shortwave
        allocate ( Cld_spec%camtsw_band           (ix, jx, kx, nswcb) )
        allocate ( Cld_spec%ncldsw_band           (ix, jx,     nswcb) )
        allocate ( Cld_spec%cld_thickness_sw_band (ix, jx, kx, nswcb) )
        allocate ( Cld_spec%iwp_sw_band           (ix, jx, kx, nswcb) )
        allocate ( Cld_spec%lwp_sw_band           (ix, jx, kx, nswcb) )
        allocate ( Cld_spec%reff_liq_sw_band      (ix, jx, kx, nswcb) )
        allocate ( Cld_spec%reff_ice_sw_band      (ix, jx, kx, nswcb) )
        do n = 1, nswcb
          Cld_spec%camtsw_band(:,:,:,n) = 0.0
          Cld_spec%ncldsw_band(:,:,n)   = 0
          Cld_spec%cld_thickness_sw_band(:,:,:,n) = 0
          Cld_spec%lwp_sw_band(:,:,:,n)      = 0.0
          Cld_spec%iwp_sw_band(:,:,:,n)      = 0.0
          Cld_spec%reff_liq_sw_band(:,:,:,n) = 10.0
          Cld_spec%reff_ice_sw_band(:,:,:,n) = 30.0
        end do

        ! longwave
        allocate ( Cld_spec%crndlw_band           (ix, jx, kx, nlwcb) )
        allocate ( Cld_spec%nrndlw_band           (ix, jx,     nlwcb) )
        allocate ( Cld_spec%cld_thickness_lw_band (ix, jx, kx, nlwcb) )
        allocate (Cld_spec%iwp_lw_band            (ix, jx, kx, nlwcb) )
        allocate (Cld_spec%lwp_lw_band            (ix, jx, kx, nlwcb) )
        allocate (Cld_spec%reff_liq_lw_band       (ix, jx, kx, nlwcb) )
        allocate (Cld_spec%reff_ice_lw_band       (ix, jx, kx, nlwcb) )
        do n = 1, nlwcb
          Cld_spec%crndlw_band(:,:,:,n) = 0.0
          Cld_spec%nrndlw_band(:,:,n)   = 0
          Cld_spec%cld_thickness_lw_band(:,:,:,n) = 0
          Cld_spec%lwp_lw_band(:,:,:,n)       = 0.0
          Cld_spec%iwp_lw_band(:,:,:,n)       = 0.0
          Cld_spec%reff_liq_lw_band(:,:,:,n)  = 10.0
          Cld_spec%reff_ice_lw_band (:,:,:,n) = 30.0
        end do

        allocate (Cld_spec%stoch_cloud_type (ix, jx, kx, nswcb + nlwcb) )
        Cld_spec%stoch_cloud_type = 0
      endif

!--------------------------------------------------------------------
!    allocate and initialize various arrays that are used by one or
!    another cloud scheme to specify the cloud locations and amounts.
!    initialization provides values consistent with the absence of
!    cloud, with the exception of the particle size fields which are
!    set to small, non-zero values.
!---------------------------------------------------------------------
      allocate (Cld_spec%hi_cloud       (ix, jx, kx) )
      allocate (Cld_spec%mid_cloud      (ix, jx, kx) )
      allocate (Cld_spec%low_cloud      (ix, jx, kx) )
      allocate (Cld_spec%ice_cloud      (ix, jx, kx) )
      allocate (Cld_spec%iwp            (ix, jx, kx) )
      allocate (Cld_spec%lwp            (ix, jx, kx) )
      allocate (Cld_spec%reff_liq       (ix, jx, kx) )
      allocate (Cld_spec%reff_ice       (ix, jx, kx) )
      allocate (Cld_spec%reff_liq_lim   (ix, jx, kx) )
      allocate (Cld_spec%reff_ice_lim   (ix, jx, kx) )
      allocate (Cld_spec%reff_liq_micro (ix, jx, kx) )
      allocate (Cld_spec%reff_ice_micro (ix, jx, kx) )
!BW   allocate (Cld_spec%tau            (ix, jx, kx, num_slingo_bands) )
      allocate (Cld_spec%liq_frac       (ix, jx, kx) )
      allocate (Cld_spec%cld_thickness  (ix, jx, kx) )
      allocate (Cld_spec%cloud_water    (ix, jx, kx) )
      allocate (Cld_spec%cloud_ice      (ix, jx, kx) )
      allocate (Cld_spec%cloud_area     (ix, jx, kx) )
      allocate (Cld_spec%cloud_droplet  (ix, jx, kx) )
      allocate (Cld_spec%cloud_ice_num  (ix, jx, kx) )
      allocate (Cld_spec%snow           (ix, jx, kx) )
      allocate (Cld_spec%rain           (ix, jx, kx) )
      allocate (Cld_spec%snow_size      (ix, jx, kx) )
      allocate (Cld_spec%rain_size      (ix, jx, kx) )

      Cld_spec%hi_cloud (:,:,:)      = .false.
      Cld_spec%mid_cloud(:,:,:)      = .false.
      Cld_spec%low_cloud(:,:,:)      = .false.
      Cld_spec%ice_cloud(:,:,:)      = .false.
      Cld_spec%lwp(:,:,:)            = 0.0
      Cld_spec%iwp(:,:,:)            = 0.0
      Cld_spec%reff_liq(:,:,:)       = 10.0
      Cld_spec%reff_ice(:,:,:)       = 30.0
      Cld_spec%reff_liq_lim(:,:,:)   = 10.0
      Cld_spec%reff_ice_lim(:,:,:)   = 30.0
      Cld_spec%reff_liq_micro(:,:,:) = 10.0
      Cld_spec%reff_ice_micro(:,:,:) = 30.0
      Cld_spec%liq_frac(:,:,:)       = 0.0
      Cld_spec%cld_thickness(:,:,:)  = 0
      Cld_spec%cloud_water(:,:,:)    = 0.0
      Cld_spec%cloud_ice(:,:,:)      = 0.0
      Cld_spec%cloud_area(:,:,:)     = 0.0
      Cld_spec%cloud_droplet(:,:,:)  = 0.0
      Cld_spec%cloud_ice_num(:,:,:)  = 0.0
      Cld_spec%snow(:,:,:)           = 0.0
      Cld_spec%rain(:,:,:)           = 0.0
      Cld_spec%snow_size(:,:,:)      = 0.0
      Cld_spec%rain_size(:,:,:)      = 0.0
!BW   do n=1, num_slingo_bands
!BW     Cld_spec%tau(:,:,:,n)  = 0.0
!BW   end do

!--------------------------------------------------------------------

end subroutine cloud_spec_alloc

!####################################################################
! <SUBROUTINE NAME="cloud_spec_dealloc">
!  <OVERVIEW>
!    cloud_spec_dealloc deallocates the component arrays of the 
!    cld_specification_type structure Cld_spec
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_spec_dealloc deallocates the component arrays of the 
!    cld_specification_type structure Cld_spec
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec_dealloc (Cld_spec, Lsc_microphys)
!  </TEMPLATE>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
! </SUBROUTINE>
! 
subroutine cloud_spec_dealloc (Cld_spec, Cldrad_control)

!---------------------------------------------------------------------
!    cloud_spec_dealloc deallocates the component arrays of the 
!    cld_specification_type structure Cld_spec
!----------------------------------------------------------------------

class(cld_specification_type), intent(inout) :: Cld_spec
type(cloudrad_control_type),   intent(in)    :: Cldrad_control

!----------------------------------------------------------------------
!    deallocate the array elements of Cld_spec.
!----------------------------------------------------------------------
      deallocate (Cld_spec%camtsw         )
      deallocate (Cld_spec%cmxolw         )
      deallocate (Cld_spec%crndlw         )
      deallocate (Cld_spec%ncldsw         )
      deallocate (Cld_spec%nmxolw         )
      deallocate (Cld_spec%nrndlw         )
      if (Cldrad_control%do_stochastic_clouds) then
        deallocate (Cld_spec%camtsw_band    )
        deallocate (Cld_spec%ncldsw_band    )
        deallocate (Cld_spec%cld_thickness_sw_band  )
        deallocate (Cld_spec%lwp_sw_band            )
        deallocate (Cld_spec%iwp_sw_band            )
        deallocate (Cld_spec%reff_liq_sw_band       )
        deallocate (Cld_spec%reff_ice_sw_band       )
        deallocate (Cld_spec%stoch_cloud_type       )
      endif
      if (Cldrad_control%do_stochastic_clouds) then
        deallocate (Cld_spec%crndlw_band    )
        deallocate (Cld_spec%nrndlw_band    )
        deallocate (Cld_spec%cld_thickness_lw_band  )
        deallocate (Cld_spec%lwp_lw_band            )
        deallocate (Cld_spec%iwp_lw_band            )
        deallocate (Cld_spec%reff_liq_lw_band       )
        deallocate (Cld_spec%reff_ice_lw_band       )
      endif
!BW   deallocate (Cld_spec%tau            )
      deallocate (Cld_spec%lwp            )
      deallocate (Cld_spec%iwp            )
      deallocate (Cld_spec%reff_liq       )
      deallocate (Cld_spec%reff_ice       )
      deallocate (Cld_spec%reff_liq_lim   )
      deallocate (Cld_spec%reff_ice_lim   )
      deallocate (Cld_spec%reff_liq_micro )
      deallocate (Cld_spec%reff_ice_micro )
      deallocate (Cld_spec%liq_frac       )
      deallocate (Cld_spec%cld_thickness  )
      deallocate (Cld_spec%hi_cloud       )
      deallocate (Cld_spec%mid_cloud      )
      deallocate (Cld_spec%low_cloud      )
      deallocate (Cld_spec%ice_cloud      )
      deallocate (Cld_spec%cloud_water    )
      deallocate (Cld_spec%cloud_ice      )
      deallocate (Cld_spec%cloud_area     )
      deallocate (Cld_spec%cloud_droplet  )
      deallocate (Cld_spec%cloud_ice_num  )
      deallocate (Cld_spec%snow           )
      deallocate (Cld_spec%rain           )
      deallocate (Cld_spec%snow_size      )
      deallocate (Cld_spec%rain_size      )

!--------------------------------------------------------------------

end subroutine cloud_spec_dealloc

!####################################################################

subroutine microphysics_alloc (Cloud_microphys, ix, jx, kx, &
                               scheme_name, Cldrad_control)

class(microphysics_type),    intent(inout) :: Cloud_microphys
integer,                     intent(in)    :: ix, jx, kx
character(len=*),            intent(in)    :: scheme_name
type(cloudrad_control_type), intent(in)    :: Cldrad_control

!--------------------------------------------------------------------
! local variables

      integer  :: nb, nswcb, nlwcb

!--------------------------------------------------------------------

      nswcb = Cldrad_control%num_sw_cloud_bands
      nlwcb = Cldrad_control%num_lw_cloud_bands

      Cloud_microphys%scheme_name = trim(scheme_name)
      Cloud_microphys%use_for_diag = .false.

!---------------------------------------------------------------------
!    allocate the arrays defining the microphysical parameters
!    common to all cloud schemes, including any precipitation fields
!    and the total cloud fraction. concentrations and fractions are
!    initialized to 0.0, and the effective sizes are set to small
!    numbers to avoid potential divides by zero.
!---------------------------------------------------------------------
      allocate (Cloud_microphys%conc_drop  (ix, jx, kx) )
      allocate (Cloud_microphys%conc_ice   (ix, jx, kx) )
      allocate (Cloud_microphys%conc_rain  (ix, jx, kx) )
      allocate (Cloud_microphys%conc_snow  (ix, jx, kx) )
      allocate (Cloud_microphys%size_drop  (ix, jx, kx) )
      allocate (Cloud_microphys%size_ice   (ix, jx, kx) )
      allocate (Cloud_microphys%size_rain  (ix, jx, kx) )
      allocate (Cloud_microphys%size_snow  (ix, jx, kx) )
      allocate (Cloud_microphys%cldamt     (ix, jx, kx) )
      allocate (Cloud_microphys%droplet_number (ix, jx, kx) )
      Cloud_microphys%conc_drop(:,:,:) = 0.
      Cloud_microphys%conc_ice(:,:,:)  = 0.
      Cloud_microphys%conc_rain(:,:,:) = 0.
      Cloud_microphys%conc_snow(:,:,:) = 0.
      Cloud_microphys%size_drop(:,:,:) = 1.0e-20
      Cloud_microphys%size_ice(:,:,:)  = 1.0e-20
      Cloud_microphys%size_rain(:,:,:) = 1.0e-20
      Cloud_microphys%size_snow(:,:,:) = 1.0e-20
      Cloud_microphys%cldamt(:,:,:)    = 0.0
      Cloud_microphys%droplet_number(:,:,:) = 0.0

!---------------------------------------------------------------------
!    allocate the arrays unique to the large-scale stratiform clouds 
!---------------------------------------------------------------------
      if (trim(scheme_name) == 'strat_cloud' .or. trim(scheme_name) == 'diag') then

        allocate (Cloud_microphys%ice_number (ix, jx, kx) )
        Cloud_microphys%ice_number(:,:,:)= 0.0

    !--------------------------------------------------
    !    allocate the arrays for stochastic clouds
    !--------------------------------------------------
        if (Cldrad_control%do_stochastic_clouds) then
          allocate (Cloud_microphys%stoch_conc_drop (ix, jx, kx, nlwcb + nswcb) )
          allocate (Cloud_microphys%stoch_conc_ice  (ix, jx, kx, nlwcb + nswcb) )
          allocate (Cloud_microphys%stoch_size_drop (ix, jx, kx, nlwcb + nswcb) )
          allocate (Cloud_microphys%stoch_size_ice  (ix, jx, kx, nlwcb + nswcb) )
          allocate (Cloud_microphys%stoch_cldamt    (ix, jx, kx, nlwcb + nswcb) )
          allocate (Cloud_microphys%stoch_droplet_number (ix, jx, kx, nlwcb + nswcb) )
          allocate (Cloud_microphys%stoch_ice_number     (ix, jx, kx, nlwcb + nswcb) )

          Cloud_microphys%lw_stoch_conc_drop => Cloud_microphys%stoch_conc_drop(:, :, :, 1:nlwcb)
          Cloud_microphys%sw_stoch_conc_drop => Cloud_microphys%stoch_conc_drop(:, :, :, nlwcb+1:)
          Cloud_microphys%lw_stoch_conc_ice  => Cloud_microphys%stoch_conc_ice (:, :, :, 1:nlwcb)
          Cloud_microphys%sw_stoch_conc_ice  => Cloud_microphys%stoch_conc_ice (:, :, :, nlwcb+1:)
          Cloud_microphys%lw_stoch_size_drop => Cloud_microphys%stoch_size_drop(:, :, :, 1:nlwcb)
          Cloud_microphys%sw_stoch_size_drop => Cloud_microphys%stoch_size_drop(:, :, :, nlwcb+1:)
          Cloud_microphys%lw_stoch_size_ice  => Cloud_microphys%stoch_size_ice (:, :, :, 1:nlwcb)
          Cloud_microphys%sw_stoch_size_ice  => Cloud_microphys%stoch_size_ice (:, :, :, nlwcb+1:)
          Cloud_microphys%lw_stoch_cldamt    => Cloud_microphys%stoch_cldamt   (:, :, :, 1:nlwcb)
          Cloud_microphys%sw_stoch_cldamt    => Cloud_microphys%stoch_cldamt   (:, :, :, nlwcb+1:)
          Cloud_microphys%lw_stoch_droplet_number    => Cloud_microphys%stoch_droplet_number   (:, :, :, 1:nlwcb)
          Cloud_microphys%sw_stoch_droplet_number    => Cloud_microphys%stoch_droplet_number   (:, :, :, nlwcb+1:)
          Cloud_microphys%lw_stoch_ice_number    => Cloud_microphys%stoch_ice_number   (:, :, :, 1:nlwcb)
          Cloud_microphys%sw_stoch_ice_number    => Cloud_microphys%stoch_ice_number   (:, :, :, nlwcb+1:)

          do nb = 1, nlwcb + nswcb
            Cloud_microphys%stoch_conc_drop(:,:,:,nb) = 0.
            Cloud_microphys%stoch_conc_ice(:,:,:,nb)  = 0.
            Cloud_microphys%stoch_size_drop(:,:,:,nb) = 1.0e-20
            Cloud_microphys%stoch_size_ice(:,:,:,nb)  = 1.0e-20
            Cloud_microphys%stoch_cldamt(:,:,:,nb)    = 0.0
            Cloud_microphys%stoch_droplet_number(:,:,:,nb) = 0.0
            Cloud_microphys%stoch_ice_number(:,:,:,nb)= 0.0
          enddo

        ! for diagnostics of all schemes
          if (trim(scheme_name) == 'diag') then
            allocate (Cloud_microphys%stoch_cloud_type (ix, jx, kx, nlwcb + nswcb) )
            Cloud_microphys%use_for_diag = .true.
          endif

        endif ! do_stochastic_clouds

      ! allocate arrays to store large-scale cloud diagnostics
        if (trim(scheme_name) == 'diag') then
          allocate(Cloud_microphys%lsc_cldamt        (ix, jx, kx)) 
          allocate(Cloud_microphys%lsc_conc_drop     (ix, jx, kx)) 
          allocate(Cloud_microphys%lsc_size_drop     (ix, jx, kx)) 
          allocate(Cloud_microphys%lsc_droplet_number(ix, jx, kx)) 
          Cloud_microphys%lsc_cldamt    = 0.0
          Cloud_microphys%lsc_conc_drop = 0.0
          Cloud_microphys%lsc_size_drop = 0.0
          Cloud_microphys%lsc_droplet_number = 0.0
        endif

      endif ! strat_cloud

!--------------------------------------------------------------------

end subroutine microphysics_alloc

!####################################################################

subroutine microphysics_alloc_diag (Model_microphys, ix, jx, kx, &
                                    scheme_names, Cldrad_control)

class(microphysics_type),       intent(inout) :: Model_microphys
integer,                        intent(in)    :: ix, jx, kx
character(len=*), dimension(:), intent(in)    :: scheme_names
type(cloudrad_control_type),    intent(in)    :: Cldrad_control
!--------------------------------------------------------------------
!  local variables

      integer :: n
      character(len=64) :: all_scheme_names

!----------------------------------------------------------------------
!    for the diagnostic microphysics_type the cloud scheme name
!    contains a comma-separated list all cloud schemes names
!    this is needed for determining stochastic cloud types in cosp
!----------------------------------------------------------------------
      all_scheme_names = ''
      if (size(scheme_names(:)) > 0) then
        all_scheme_names = scheme_names(1)
        do n = 2, size(scheme_names(:))
          all_scheme_names = trim(all_scheme_names)//','//trim(scheme_names(n))
        enddo
      endif

      call Model_microphys%alloc( ix, jx, kx, 'diag', Cldrad_control)

      Model_microphys%scheme_name = all_scheme_names

!--------------------------------------------------------------------

end subroutine microphysics_alloc_diag

!####################################################################

subroutine microphysics_dealloc (Cloud_microphys, Cldrad_control)

class(microphysics_type),    intent(inout) :: Cloud_microphys
type(cloudrad_control_type), intent(in)    :: Cldrad_control
!--------------------------------------------------------------------
!    deallocates the component arrays of the microphysics_type
!    data structure
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    deallocate the elements common to all cloud types
!--------------------------------------------------------------------
      deallocate (Cloud_microphys%conc_drop   )
      deallocate (Cloud_microphys%conc_ice    )
      deallocate (Cloud_microphys%conc_rain   )
      deallocate (Cloud_microphys%conc_snow   )
      deallocate (Cloud_microphys%size_drop   )
      deallocate (Cloud_microphys%size_ice    )
      deallocate (Cloud_microphys%size_rain   )
      deallocate (Cloud_microphys%size_snow   )
      deallocate (Cloud_microphys%cldamt      )
      deallocate (Cloud_microphys%droplet_number )

!---------------------------------------------------------------------
!    deallocate the arrays unique to the large-scale stratiform clouds 
!---------------------------------------------------------------------
      if (trim(Cloud_microphys%scheme_name) == 'strat_cloud' .or. Cloud_microphys%use_for_diag) then

        deallocate(Cloud_microphys%ice_number)

    !--------------------------------------------------
    !    deallocate the arrays for stochastic clouds
    !--------------------------------------------------
        if (Cldrad_control%do_stochastic_clouds) then
          deallocate (Cloud_microphys%stoch_conc_drop)
          deallocate (Cloud_microphys%stoch_conc_ice)
          deallocate (Cloud_microphys%stoch_size_drop)
          deallocate (Cloud_microphys%stoch_size_ice)
          deallocate (Cloud_microphys%stoch_cldamt)
          deallocate (Cloud_microphys%stoch_droplet_number)
          deallocate (Cloud_microphys%stoch_ice_number)

          nullify (Cloud_microphys%lw_stoch_conc_drop   )
          nullify (Cloud_microphys%lw_stoch_conc_ice    )
          nullify (Cloud_microphys%lw_stoch_size_drop   )
          nullify (Cloud_microphys%lw_stoch_size_ice    )
          nullify (Cloud_microphys%lw_stoch_cldamt      )
          nullify (Cloud_microphys%lw_stoch_droplet_number)
          nullify (Cloud_microphys%lw_stoch_ice_number  )

          nullify (Cloud_microphys%sw_stoch_conc_drop   )
          nullify (Cloud_microphys%sw_stoch_conc_ice    )
          nullify (Cloud_microphys%sw_stoch_size_drop   )
          nullify (Cloud_microphys%sw_stoch_size_ice    )
          nullify (Cloud_microphys%sw_stoch_cldamt      )
          nullify (Cloud_microphys%sw_stoch_droplet_number)
          nullify (Cloud_microphys%sw_stoch_ice_number  )

          if (Cloud_microphys%use_for_diag) then
            deallocate (Cloud_microphys%stoch_cloud_type)
          endif
        endif ! do_stochastic_clouds

      ! deallocate arrays to store large-scale cloud diagnostics
        if (Cloud_microphys%use_for_diag) then
          deallocate(Cloud_microphys%lsc_cldamt        )
          deallocate(Cloud_microphys%lsc_conc_drop     )
          deallocate(Cloud_microphys%lsc_size_drop     )
          deallocate(Cloud_microphys%lsc_droplet_number)
        endif

      endif ! strat_cloud

      Cloud_microphys%scheme_name = ' '
      Cloud_microphys%use_for_diag = .false.

!--------------------------------------------------------------------

end subroutine microphysics_dealloc

!####################################################################

subroutine cldrad_properties_alloc (Cldrad_props, ix, jx, kx, &
                                    Cldrad_control)

class(cldrad_properties_type),   intent(inout) :: Cldrad_props
integer,                         intent(in)    :: ix, jx, kx
type(cloudrad_control_type),     intent(in)    :: Cldrad_control

integer :: nlwcb, nswcb
!--------------------------------------------------------------------

      nlwcb = Cldrad_control%num_lw_cloud_bands
      nswcb = Cldrad_control%num_sw_cloud_bands

!-------------------------------------------------------------------
!    allocate the arrays used to define the longwave cloud radiative
!    properties. initialize to appropriate non-cloudy values.
!--------------------------------------------------------------------
      if (Cldrad_control%do_ica_calcs) then
        allocate (Cldrad_props%emmxolw  (ix, jx, kx, nlwcb, nlwcb) )
        allocate (Cldrad_props%emrndlw  (ix, jx, kx, nlwcb, nlwcb) )
        allocate (Cldrad_props%abscoeff (ix, jx, kx, nlwcb, nlwcb) )
        allocate (Cldrad_props%cldemiss (ix, jx, kx, nlwcb, nlwcb) )
      else
        allocate (Cldrad_props%emmxolw  (ix, jx, kx, nlwcb, 1) )
        allocate (Cldrad_props%emrndlw  (ix, jx, kx, nlwcb, 1) )
        allocate (Cldrad_props%abscoeff (ix, jx, kx, nlwcb, 1) )
        allocate (Cldrad_props%cldemiss (ix, jx, kx, nlwcb, 1) )
      endif
      Cldrad_props%emmxolw  = 1.0
      Cldrad_props%emrndlw  = 1.0
      Cldrad_props%abscoeff = 0.0
      Cldrad_props%cldemiss = 0.0

!---------------------------------------------------------------------
!    allocate and initialize the microphysically-based shortwave cloud 
!    radiative properties.
!---------------------------------------------------------------------
      if (Cldrad_control%do_ica_calcs) then
        allocate (Cldrad_props%cldext  (ix, jx, kx, nswcb, nswcb) )
        allocate (Cldrad_props%cldsct  (ix, jx, kx, nswcb, nswcb) )
        allocate (Cldrad_props%cldasymm(ix, jx, kx, nswcb, nswcb) )
      else
        allocate (Cldrad_props%cldext  (ix, jx, kx, nswcb, 1))
        allocate (Cldrad_props%cldsct  (ix, jx, kx, nswcb, 1))
        allocate (Cldrad_props%cldasymm(ix, jx, kx, nswcb, 1))
      endif
      Cldrad_props%cldsct   = 0.0
      Cldrad_props%cldext   = 0.0
      Cldrad_props%cldasymm = 1.0

!--------------------------------------------------------------------

end subroutine cldrad_properties_alloc

!####################################################################

subroutine cldrad_properties_dealloc (Cldrad_props)

class(cldrad_properties_type),   intent(inout) :: Cldrad_props
!--------------------------------------------------------------------

     ! longwave
      deallocate (Cldrad_props%emmxolw   )
      deallocate (Cldrad_props%emrndlw   )
      deallocate (Cldrad_props%abscoeff  )
      deallocate (Cldrad_props%cldemiss  )
     ! shortwave
      deallocate (Cldrad_props%cldext    )
      deallocate (Cldrad_props%cldasymm  )
      deallocate (Cldrad_props%cldsct    )

!--------------------------------------------------------------------

end subroutine cldrad_properties_dealloc

!####################################################################

subroutine microrad_properties_alloc (Microrad_props, ix, jx, kx, &
                                      Cldrad_control)

type(microrad_properties_type), intent(inout) :: Microrad_props
integer,                         intent(in)    :: ix, jx, kx
type(cloudrad_control_type),     intent(in)    :: Cldrad_control
!--------------------------------------------------------------------

      allocate (Microrad_props%cldext   (ix, jx, kx, Cldrad_control%num_sw_cloud_bands))
      allocate (Microrad_props%cldsct   (ix, jx, kx, Cldrad_control%num_sw_cloud_bands))
      allocate (Microrad_props%cldasymm (ix, jx, kx, Cldrad_control%num_sw_cloud_bands))
      allocate (Microrad_props%abscoeff (ix, jx, kx, Cldrad_control%num_lw_cloud_bands))

!--------------------------------------------------------------------

end subroutine microrad_properties_alloc

!####################################################################

subroutine microrad_properties_dealloc (Microrad_props)

class(microrad_properties_type), intent(inout) :: Microrad_props
!--------------------------------------------------------------------

      deallocate (Microrad_props%cldext  )
      deallocate (Microrad_props%cldsct  )
      deallocate (Microrad_props%cldasymm)
      deallocate (Microrad_props%abscoeff)

!--------------------------------------------------------------------

end subroutine microrad_properties_dealloc

!####################################################################

end module cloudrad_types_mod

