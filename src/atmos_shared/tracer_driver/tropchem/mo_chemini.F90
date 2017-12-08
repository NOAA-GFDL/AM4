
      module mo_chemini_mod

implicit none
      private
      public :: chemini

character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
logical                       :: module_is_initialized = .false.

      contains

      subroutine chemini( file_jval_lut, file_jval_lut_min, use_tdep_jvals, &
                          o3_column_top, jno_scale_factor, verbose, &
                          retain_cm3_bugs, do_fastjx_photo , trop_option)
!-----------------------------------------------------------------------
!       ... Chemistry module intialization
!-----------------------------------------------------------------------

      use MO_PHOTO_MOD,      only : prate_init
      use mo_chem_utls_mod,  only : chem_utls_init
      use mo_usrrxt_mod,     only : usrrxt_init
#ifndef AM3_CHEM
      use CHEM_MODS_mod,     only : grpcnt, clscnt1, clscnt4, clscnt5, chem_mods_init
#else
      use AM3_CHEM_MODS_mod, only : grpcnt, clscnt1, clscnt4, clscnt5, chem_mods_init
#endif
      use MO_EXP_SOL_mod,    only : exp_slv_init
      use MO_IMP_SOL_mod,    only : imp_slv_init
      use MO_RODAS_SOL_mod,  only : rodas_slv_init

      use MO_READ_SIM_CHM_mod, only : read_sim_chm
#ifndef AM3_CHEM
      use mo_fphoto_mod,     only : fprate_init
#else
      use AM3_fphoto_mod,     only : fprate_init
#endif
      use tropchem_types_mod, only : tropchem_opt

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: file_jval_lut, &
                                      file_jval_lut_min
      logical,          intent(in) :: use_tdep_jvals
      real,             intent(in) :: o3_column_top, &
                                      jno_scale_factor
      integer,          intent(in) :: verbose
      logical,          intent(in) :: retain_cm3_bugs
      logical,          intent(in) :: do_fastjx_photo
      type(tropchem_opt), intent(in) :: trop_option        

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      character(len=80) ::   lpath
      character(len=80) ::   mspath
      character(len=32) ::   filename, filename_solarmin
      
      character(len=128) ::  sim_data_flsp
      integer :: sim_file_cnt

!-----------------------------------------------------------------------
!       ... Allocate variables
!-----------------------------------------------------------------------
      call chem_mods_init

!-----------------------------------------------------------------------
!       ... Read sim.dat
!-----------------------------------------------------------------------
!     sim_data_flsp = 'INPUT/sim.dat'
      sim_data_flsp = 'INPUT/' // TRIM(trop_option%sim_data_flsp)
      call read_sim_chm( sim_data_flsp, sim_file_cnt )

!-----------------------------------------------------------------------
!       ... Diagnostics initialization
!-----------------------------------------------------------------------
!     call diags_init( tracnam, plonl, platl, pplon )

!-----------------------------------------------------------------------
!     ... initialize fast-jx photo
!-----------------------------------------------------------------------    
      if (do_fastjx_photo) then
         call fprate_init (o3_column_top)
      else

!-----------------------------------------------------------------------
!       ... Initialize photorate module
!-----------------------------------------------------------------------
!     filename = photo_flsp%nl_filename
!     lpath    = photo_flsp%local_path
!     mspath   = photo_flsp%remote_path
         lpath = 'INPUT/'
         filename = TRIM(file_jval_lut)
         filename_solarmin = TRIM(file_jval_lut_min)
         call prate_init( filename, filename_solarmin, lpath, mspath, &
                          use_tdep_jvals, o3_column_top, jno_scale_factor, &
                          retain_cm3_bugs )
      end if

!-----------------------------------------------------------------------
!       ... Read time-independent airplane emissions
!-----------------------------------------------------------------------
!     emires = emis_flsp%hor_res
!     if( emires(1:1) /= '.' ) then
!        emires = '.' // emires
!     end if
!     lpath    = emis_flsp%local_path
!     mspath   = emis_flsp%remote_path
!     filename = 'emissions.aircraft' // TRIM(emires) // '.nc'
!     call airpl_src( filename, lpath, mspath, plonl, platl, pplon )

!-----------------------------------------------------------------------
!       ... Initialize the chem utils module
!-----------------------------------------------------------------------
      call chem_utls_init (retain_cm3_bugs)

!-----------------------------------------------------------------------
!       ... Read time-dependent surface flux dataset
!-----------------------------------------------------------------------
!     call srf_emis_init( plonl, platl, pplon )

!-----------------------------------------------------------------------
!       ... Intialize the het rates module
!-----------------------------------------------------------------------
!     call sethet_init

!-----------------------------------------------------------------------
!       ... Intialize the ext frcing module
!-----------------------------------------------------------------------
!     call setext_init

!-----------------------------------------------------------------------
!       ... Intialize the rxt rate constant module
!-----------------------------------------------------------------------
      call usrrxt_init( verbose , trop_option)

!-----------------------------------------------------------------------
!       ... Intialize the grp ratios module
!-----------------------------------------------------------------------
!     call set_grp_ratios_init

!-----------------------------------------------------------------------
!       ... Read time-dependent surface variables dataset
!-----------------------------------------------------------------------
!     surfres = surf_flsp%hor_res
!     if( surfres(1:1) /= '.' ) then
!        surfres = '.' // surfres
!     end if
!     filename = 'surfdata' // TRIM(surfres) // '.nc'
!     lpath    = surf_flsp%local_path
!     mspath   = surf_flsp%remote_path
!     call surf_init( filename, lpath, mspath, plonl, platl, pplon )

!-----------------------------------------------------------------------
!       ... Read time-dependent upper boundary values
!-----------------------------------------------------------------------
!     filename = ubc_flsp%nl_filename
!     lpath    = ubc_flsp%local_path
!     mspath   = ubc_flsp%remote_path
!     call ub_init( platl, filename, lpath, mspath )

!-----------------------------------------------------------------------
!       ... Read time-dependent sulfate dataset
!           NOTE : This is now a netcdf dataset
!-----------------------------------------------------------------------
!     filename = 'sulfate.M1.nc'
!     lpath    = sulf_flsp%local_path
!     mspath   = sulf_flsp%remote_path
!     call sulf_init( plonl, platl, pplon, filename, lpath, mspath )

      if( clscnt1 > 0 ) then
!-----------------------------------------------------------------------
!       ... Initialize the explicit solver
!-----------------------------------------------------------------------
         call exp_slv_init
      end if
      if( clscnt4 > 0 ) then
!-----------------------------------------------------------------------
!       ... Initialize the implicit solver
!-----------------------------------------------------------------------
         call imp_slv_init( verbose, retain_cm3_bugs )
      end if
      if( clscnt5 > 0 ) then
!-----------------------------------------------------------------------
!       ... Initialize the implicit solver
!-----------------------------------------------------------------------
         call rodas_slv_init (retain_cm3_bugs)
      end if


      end subroutine chemini

      end module mo_chemini_mod
