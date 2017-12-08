!FDOC_TAG_GFDL
                        module lscloud_debug_mod


use mpp_mod,               only : input_nml_file
use fms_mod,               only : file_exist, open_namelist_file,&
                                  open_file, error_mesg, FATAL, NOTE, &
                                  mpp_pe, mpp_root_pe, close_file, &
                                  stdlog, check_nml_error,  &
                                  write_version_number, stdout
use lscloud_types_mod,     only : lscloud_types_init, lscloud_debug_type, &
                                  cloud_processes_type, atmos_state_type, &
                                  cloud_state_type, precip_state_type
use moist_proc_utils_mod,  only : mp_input_type, mp_conv2ls_type, &
                                  mp_removal_type
use check_nan_mod,         only : check_nan_init, check_nan


implicit none
private

!------------------------------------------------------------------------
!---interfaces-----------------------------------------------------------

public                                          &
          lscloud_debug_init, lscloud_debug_setup, lscloud_debug, &
          macrophysics_debug, microphysics_debug1, microphysics_debug2,  &
          aerosol_cloud_debug1, aerosol_cloud_debug2, &
          write_debug_output, record_micro_refusal, output_refusals


!------------------------------------------------------------------------
!---version number-------------------------------------------------------

Character(len=128) :: Version = '$Id: $'
Character(len=128) :: Tagname = '$Name: $'

!------------------------------------------------------------------------
!---namelist-------------------------------------------------------------

!  <DATA NAME="num_strat_pts" TYPE="integer"DEFAULT="0">
!   number of grid points where instantaneous output will be saved to
!   file strat.data.
!   num_strat_pts must be <= max_strat_pts.
!   NOTE: this option is currently not supported.
!  </DATA>
!  <DATA NAME="strat_pts" TYPE="integer" DEFAULT="0,0">
!   "num_strat_pts" pairs of grid indices, i.e., the global indices for i,j.
!   NOTE: this option is currently not supported.
!  </DATA>
!  <DATA NAME="debugo" TYPE="logical" DEFAULT=".false.">
!   should we activated debug diagnostics ?
!  </DATA>
!  <DATA NAME="debugo0" TYPE="logical" DEFAULT=".false.">
!   should we activated debug diagnostics  -- not fully implemented?
!  </DATA>
!  <DATA NAME="debugo4" TYPE="logical" DEFAULT=".false.">
!   if true, nrefuse is output
!  </DATA>
!  <DATA NAME="isamp" TYPE="integer" DEFAULT="1">
!   i coordinate of grid point for which debug data is to be written
!  </DATA>
!  <DATA NAME="jsamp" TYPE="integer" DEFAULT="1">
!   j coordinate of grid point for which debug data is to be written
!  </DATA>
!  <DATA NAME="ksamp" TYPE="integer" DEFAULT="1">
!   k coordinate of grid point for which debug data is to be written
!  </DATA>

integer, PARAMETER :: max_strat_pts = 5
integer                             :: num_strat_pts  =  0
integer,dimension(2,max_strat_pts)  :: strat_pts = 0
logical                             :: debugo = .false.
logical                             :: debugo0 = .false.
logical                             :: debugo4 = .false.
integer                             :: isamp = 1
integer                             :: jsamp = 1
integer                             :: ksamp = 1

namelist / lscloud_debug_nml /        &
              num_strat_pts, strat_pts, &
              debugo, debugo0, debugo4, isamp, jsamp, ksamp


type(lscloud_debug_type), save   :: Debug
character(len=64)                :: otname =  'debug_strt_cld.txt'

logical   :: module_is_initialized = .false.




                             CONTAINS



!#######################################################################


subroutine lscloud_debug_init (debug_out) 

logical, intent(out) :: debug_out

!------------------------------------------------------------------------
!    initialize the lscloud_debug module.
!------------------------------------------------------------------------

      integer :: io, ierr, logunit, unit

!------------------------------------------------------------------------ 
!    if module is already initialized, return.
!------------------------------------------------------------------------ 
!     if (module_is_initialized) return
      if (module_is_initialized) then
        debug_out = Debug%debugo
        return
      endif

!----------------------------------------------------------------------- 
!    process namelist.
!-----------------------------------------------------------------------  
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=lscloud_debug_nml, iostat=io)
      ierr = check_nml_error(io,'lscloud_debug_nml')
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=lscloud_debug_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'lscloud_debug_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!----------------------------------------------------------------------- 
!    write version and namelist to stdlog.
!-----------------------------------------------------------------------  
      call write_version_number(Version, Tagname)
      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() )  &
                    write (logunit, nml=lscloud_debug_nml)

!-------------------------------------------------------------------------
!    assure consistency across logical controls. debugo must be .true. for
!    either of debugo0 or debugo4 to be .true.
!      debugo: activates this debugging module
!      debugo0: if true, limits debug output to pe0.
!      debugo4: outputs microphysics refusal stats (relevant for "mg" 
!               microphysics only)
!-------------------------------------------------------------------------
      if (debugo0 .and. (.not. debugo))  &
        call error_mesg ('ls_cloud_debug_init', &
          'must have debugo be .true. for debugo0 to be .true.', FATAL)
      if (debugo4 .and. (.not. debugo))  &
        call error_mesg ('ls_cloud_debug_init', &
          'must have debugo be .true. for debugo4 to be .true.', FATAL)
     
!------------------------------------------------------------------------
!    set up debug information for use with strat_cloud.
!------------------------------------------------------------------------
      Debug%ncall = 1
      Debug%debugo =  debugo
      Debug%debugo0 = debugo0
      Debug%debugo4 = debugo4

      if (num_strat_pts > 0) then
        allocate (Debug%strat_pts(2, num_strat_pts))
        Debug%strat_pts = strat_pts
      endif

!-------------------------------------------------------------------------
!    if debugo0 is .true, then debugo is set to .true. for only the root pe.
!-------------------------------------------------------------------------
      if ( mpp_pe() .NE. mpp_root_pe() .and. Debug%debugo0 ) then
        Debug%debugo  = .FALSE.
      endif

!----------------------------------------------------------------------
!    if debugging is active, make sure the check_nan and lscloud_types 
!    modules are initialized.
!-----------------------------------------------------------------------
      if (debugo) then
        call check_nan_init
        call lscloud_types_init
      endif

!-------------------------------------------------------------------------
!    if any debugging is desired, open a file to hold the output. 
!-------------------------------------------------------------------------
      if (Debug%debugo ) then   
        Debug%otun = open_file (otname, threading = 'multi',  &
                                                        action = 'append')
      else
        Debug%otun = 0
      endif

!-----------------------------------------------------------------------
!    define output variable indicating if debug is to be active.
!-----------------------------------------------------------------------
      debug_out = Debug%debugo
 
      module_is_initialized = .true.

!-------------------------------------------------------------------------

end subroutine lscloud_debug_init


!#########################################################################

subroutine lscloud_debug_setup (is, ie, js, je, Input_mp, Atmos_state)

!------------------------------------------------------------------------
!    subroutine lscloud_debug_setup  performs some initialization  and
!    validation of some debug variables, and outputs some fields to the
!    debug file.
!------------------------------------------------------------------------

integer,                 intent(in)    :: is, ie, js, je
type(mp_input_type),     intent(in)    :: Input_mp
type(atmos_state_type),  intent(in)    :: Atmos_state

!-----------------------------------------------------------------------
!   local variables:

      integer :: ix, jx, kx

!----------------------------------------------------------------------
!    if debug is active, perform some initialization / validation of 
!    debug parameters, and then output some requested fields.
!----------------------------------------------------------------------
      if (Debug%debugo) then

!-------------------------------------------------------------------------
!    reinitialize the counter for columns containing negative water that
!    do not engage in microphysics.
!-------------------------------------------------------------------------
        Debug%nrefuse = 0

!-------------------------------------------------------------------------
!    on the first call, be sure the values input via namelist for isamp, 
!    jsamp and ksamp are valid for the current domain decomposition; modify
!    them if needed. 
!-------------------------------------------------------------------------
        if ( Debug%ncall .eq. 1 ) then
          ix = ie - is +1
          jx = je - js +1
          kx = size(Input_mp%tin,3)

          Debug%isamp = min (ix, isamp)
          Debug%jsamp = min (jx, jsamp)
          Debug%ksamp = min (kx, ksamp)

          write(Debug%otun,*) "MODIF TEST"
          write(Debug%otun,*) "is,ie,js,je ", is, ie, js, je
          write(Debug%otun,*) "  isamp, jsamp, ksamp ",    &
                                    Debug%isamp, Debug%jsamp, Debug%ksamp
        endif

!-----------------------------------------------------------------------
!    output header with time step number to debug file, and then increment
!    the counter.
!-----------------------------------------------------------------------
        write(Debug%otun,*)"-----------------------------------------------"
        write(Debug%otun,*) "ncall ", Debug%ncall

        Debug%ncall = Debug%ncall +1

!-----------------------------------------------------------------------
!    output some requested variables at the debug point.
!-----------------------------------------------------------------------
        write (Debug%otun, *) " T, pfull ",    &
                Input_mp%tin(Debug%isamp,Debug%jsamp,Debug%ksamp),   &
                     Input_mp%pfull(Debug%isamp,Debug%jsamp,Debug%ksamp)
        write (Debug%otun, *) " aaa dqsdT, qs ",   &
                  Atmos_state%dqsdT(Debug%isamp,Debug%jsamp,Debug%ksamp),  &
                      Atmos_state%qs(Debug%isamp,Debug%jsamp,Debug%ksamp)
      endif  ! (debugo)

!----------------------------------------------------------------------

end subroutine lscloud_debug_setup



!#######################################################################

subroutine lscloud_debug (ST_out, SQ_out, Cloud_state,  Removal_mp, &
                            Precip_state, Input_mp, C2ls_mp, Atmos_State)

!------------------------------------------------------------------------
!    subroutine lscloud_debug outputs the standard debug variable suite,
!    including extreme values found in relevant lscloud output arrays, 
!    checks for NaNs within these fields, and flags any unrealistic
!    extreme values for further investiagtaion.
!------------------------------------------------------------------------

real, dimension(:,:,:),    intent(in) :: ST_out, SQ_out
type(cloud_state_type),    intent(in) :: Cloud_state
type(atmos_state_type),    intent(in) :: Atmos_state
type(mp_removal_type),     intent(in) :: Removal_mp
type(mp_input_type),       intent(in) :: Input_mp    
type(mp_conv2ls_type),     intent(in) :: C2ls_mp    
type(precip_state_type),   intent(in) :: Precip_state

!---------------------------------------------------------------------
!   local variables:

      integer, dimension(3) :: maxl, minl

!-----------------------------------------------------------------------
!    execute the subroutine if debug is active.
!-----------------------------------------------------------------------
      if (Debug%debugo) then

!------------------------------------------------------------------------
!    write the max and min values and their locations, along with 
!    information on any NaNs in the array to the debug output file.
!------------------------------------------------------------------------

!   array ST_out:
        write(Debug%otun, *)  "eee max, min ST ",   &
                                       MAXVAL(ST_out), MINVAL(ST_out)
        write(Debug%otun, *)  "eee maxloc, minloc  ",   &
                                       MAXLOC(ST_out), MINLOC(ST_out)
        call check_nan (ST_out,'ST_out')

!   array SQ_out:
        write(Debug%otun, *)  "eee max, min SQ ",   &
                                        MAXVAL(SQ_out), MINVAL(SQ_out)
        write(Debug%otun, *)  "eee maxloc, minloc  ",    &
                                       MAXLOC(SQ_out),  MINLOC(SQ_out)
        call check_nan (SQ_out,'SQ_out')

!   array SL_out:
        write(Debug%otun, *)  "eee max, min SL ",   &
                    MAXVAL(Cloud_state%SL_out), MINVAL(Cloud_state%SL_out)
        write(Debug%otun, *)  "eee maxloc, minloc  ",   &
                   MAXLOC(Cloud_state%SL_out), MINLOC(Cloud_state%SL_out)
        call check_nan   (Cloud_state%SL_out,'SL_out')

!   array SI_out:
        write(Debug%otun, *)  "eee max, min SI ",   &
                     MAXVAL(Cloud_state%SI_out), MINVAL(Cloud_state%SI_out)
        write(Debug%otun, *)  "eee maxloc, minloc  ",   &
                      MAXLOC(Cloud_state%SI_out), MINLOC(Cloud_state%SI_out)
        call check_nan   (Cloud_state%SI_out,'SI_out')

!   array SA_out:
        write(Debug%otun, *)  "eee max, min SA ",   &
                     MAXVAL(Cloud_state%SA_out), MINVAL(Cloud_state%SA_out)
        write(Debug%otun, *)  "eee maxloc, minloc  ",   &
                     MAXLOC(Cloud_state%SA_out), MINLOC(Cloud_state%SA_out)
        call check_nan   (Cloud_state%SA_out,'SA_out')

!   array SN_out:
        write(Debug%otun, *)  "eee max, min SN ",   &
                   MAXVAL(Cloud_state%SN_out),  MINVAL(Cloud_state%SN_out)
        write(Debug%otun, *)  "eee maxloc, minloc  ", &
                     MAXLOC(Cloud_state%SN_out), MINLOC(Cloud_state%SN_out)
        call check_nan   (Cloud_state%SN_out,'SN_out')

!   array SNi_out:
        write(Debug%otun, *)  "eee max, min SNi ",   &
                    MAXVAL(Cloud_state%SNi_out), MINVAL(Cloud_state%SNi_out)
        write(Debug%otun, *)  "eee maxloc, minloc  ", &
                   MAXLOC(Cloud_state%SNi_out), MINLOC(Cloud_state%SNi_out)
        call check_nan   (Cloud_state%SNi_out,'SNi_out')

!   array T + ST:
        write(Debug%otun, *)  "--"
        write(Debug%otun, *)  "eee max, min T+ST ",  &
              MAXVAL(Input_mp%tin + ST_out),  MINVAL(Input_mp%tin + ST_out)

!   array qv + SQ:
        write(Debug%otun, *)  "eee max, min qv+SQ ", &
             MAXVAL(Input_mp%qin + SQ_out), MINVAL(Input_mp%qin + SQ_out)

!   array ql + SL:
        write(Debug%otun, *)  "eee max, min ql+ SL ",   &
                     MAXVAL(Cloud_state%ql_in+Cloud_state%SL_out),  &
                               MINVAL(Cloud_state%ql_in+Cloud_state%SL_out)

!   array qi + SI:
        write(Debug%otun, *)  "eee max, min qi +SI ",  &
                      MAXVAL(Cloud_state%qi_in+ Cloud_state%SI_out), &
                             MINVAL(Cloud_state%qi_in + Cloud_state%SI_out)

!   array qa + SA:
        write(Debug%otun, *)  "eee max, min qa + SA ",  &
                           MAXVAL(Cloud_state%qa_in+Cloud_state%SA_out),  &
                              MINVAL(Cloud_state%qa_in+Cloud_state%SA_out)

!   array qn + SN:
        write(Debug%otun, *)  "eee max, min qn + SN ",  &
                         MAXVAL(Cloud_state%qn_in + Cloud_state%SN_out), &
                            MINVAL(Cloud_state%qn_in + Cloud_state%SN_out)

!   array qni + SNi:
        write(Debug%otun, *)  "eee max, min qni SNi ",   &
                        MAXVAL(Cloud_state%qni_in+ Cloud_state%SNi_out), &
                            MINVAL(Cloud_state%qni_in+Cloud_state%SNi_out)

        write(Debug%otun, *)  "--"
        write(Debug%otun, *)  "--"

!   3D rain field:
        write(Debug%otun, *)  "eee max, min rain3d ",  &
                  MAXVAL(Removal_mp%rain3d), MINVAL(Removal_mp%rain3d)
        write(debug%otun, *)  "eee maxloc, minloc  ",  &
                  MAXLOC(Removal_mp%rain3d), MINLOC(Removal_mp%rain3d)
        call check_nan   (Removal_mp%rain3d,'rain3d')

!   3D snow field:
        write(Debug%otun, *)  "eee max, min snow3d ",   &
                  MAXVAL(Removal_mp%snow3d), MINVAL(Removal_mp%snow3d)
        write(Debug%otun, *)  "eee maxloc, minloc  ",  &
                   MAXLOC(Removal_mp%snow3d), MINLOC(Removal_mp%snow3d)
        call check_nan   (Removal_mp%snow3d,'snow3d')

        write(Debug%otun, *)  "--"

!   surface precip:
        write(Debug%otun, *)  "eee max, min surfrain ",   &
               MAXVAL(Precip_state%surfrain), MINVAL(Precip_state%surfrain)
        write(Debug%otun, *)  "eee maxloc, minloc  ",    &
               MAXLOC(Precip_state%surfrain), MINLOC(Precip_state%surfrain)

!   surface snow:
        write(Debug%otun, *)  "eee max, min surfsnow ",   &
               MAXVAL(Precip_state%surfsnow), MINVAL(Precip_state%surfsnow)
        write(Debug%otun, *)  "eee maxloc, minloc  ",   &
               MAXLOC(Precip_state%surfsnow), MINLOC(Precip_state%surfsnow)
        write(Debug%otun, *)  "--"

!-----------------------------------------------------------------------
!    determine any extremes that are non-realizable.
!-----------------------------------------------------------------------

!   output qv field:
        IF ( MAXVAL(SQ_out + Input_mp%qin) .GT. 1.e-1 )  &
                                          write(debug%otun, *) " MMMMM Q1 "
        IF ( MAXVAL(SQ_out + Input_mp%qin) .LT. 0. )   &
                                          write(debug%otun, *) " MMMMM Q2 "

!   output qi field:
        IF ( MAXVAL(Cloud_state%SI_out + Cloud_state%qi_in) .LT. 0. )  &
                                         write(Debug%otun, *) " MMMMM I11 "

!   output ql field:
        IF ( MAXVAL(Cloud_state%SL_out + Cloud_state%ql_in) .LT. 0. )   &
                                        write(Debug%otun, *) " MMMMM L11 "

!   output qa field:
        IF ( MAXVAL(  &
            Cloud_state%qa_in + Cloud_state%SA_out +   &
               C2ls_mp%convective_humidity_area)  .GT. 1.00000000001) THEN
          write(Debug%otun, *) " MMMMMA1 ahuco "
          maxl =  maxloc (Cloud_state%qa_in + Cloud_state%SA_out +  &
                                          C2ls_mp%convective_humidity_area)
          write(Debug%otun, *) "  maxloc(qa+SA) ",  &
             maxloc (Cloud_state%qa_in + Cloud_state%SA_out +  &
                                         C2ls_mp%convective_humidity_area)
          write(Debug%otun, *) " qa+SA+ahuco ",    &
            Cloud_state%qa_in(maxl(1),maxl(2),maxl(3)) +     &
                Cloud_state%SA_out(maxl(1),maxl(2),maxl(3)) +&
                   C2ls_mp%convective_humidity_area(maxl(1),maxl(2),maxl(3))
          write(Debug%otun, *) " qa, ahuco, SA ",   &
             Cloud_state%qa_in(maxl(1),maxl(2),maxl(3)),   &
               C2ls_mp%convective_humidity_area(maxl(1),maxl(2),maxl(3)),  &
                           Cloud_state%SA_out(maxl(1),maxl(2),maxl(3)) 
        END IF

        IF ( MINVAL(Cloud_state%qa_in+Cloud_state%SA_out) .LT. 0.  ) THEN
          minl =  minloc(Cloud_state%qa_in+Cloud_state%SA_out)
          write(Debug%otun, *) " MMMMMA2"
          write(Debug%otun, *) "  minloc(qa+SA) ",   &
                              minloc(Cloud_state%qa_in+Cloud_state%SA_out)
          write(Debug%otun, *) " qa+SA ", &
                         Cloud_state%qa_in(minl(1),minl(2),minl(3)) +  &
                                Cloud_state%SA_out(minl(1),minl(2),minl(3))
          write(Debug%otun, *) "    qa, SA ",  &
                           Cloud_state%qa_in(minl(1),minl(2),minl(3)),  &
                                Cloud_state%SA_out(minl(1),minl(2),minl(3))
        END IF

!   output qn field:
        IF ( MINVAL(Cloud_state%qn_in+Cloud_state%SN_out) .LT. 0.  ) THEN
          write(Debug%otun, *) " MMMMMN1"
          minl =  minloc(Cloud_state%qn_in+Cloud_state%SN_out)
          write(Debug%otun, *) "  minloc(qn+SN) ",  &
                           minloc(Cloud_state%qn_in+Cloud_state%SN_out)
          write(Debug%otun, *) "    qn, SN ",   &
                          Cloud_state%qn_in(minl(1),minl(2),minl(3)),  &
                          Cloud_state%SN_out(minl(1),minl(2),minl(3)) 
        END IF

!   output qi field:
        IF ( MAXVAL(Cloud_state%qi_in +    &
                                   Cloud_state%SI_out) .GT. 1.e-2 )   then
          maxl = MAXLOC(Cloud_state%qi_in + Cloud_state%SI_out)
          write(Debug%otun, *) " MMMMMII"
          write(debug%otun, *) "  maxloc(qi+SI) ",    &
                              maxloc(Cloud_state%qi_in+Cloud_state%SI_out)
          write(Debug%otun, *) "  qi+SI " ,    &
                           Cloud_state%qi_in(maxl(1),maxl(2),maxl(3))+  &
                                Cloud_state%SI_out(maxl(1),maxl(2),maxl(3))
          write(Debug%otun, *) "  qi, SI " ,    &
                      Cloud_state%qi_in(maxl(1),maxl(2),maxl(3)),  &
                               Cloud_state%SI_out(maxl(1),maxl(2),maxl(3))
          write(Debug%otun, *) "  T ",    &
                     Input_mp%tin(maxl(1),maxl(2),maxl(3)),   &
                               Cloud_state%SI_out(maxl(1),maxl(2),maxl(3))
        END IF

!   output ql field:
        IF ( MAXVAL(Cloud_state%ql_in + Cloud_state%SL_out) .GT. 1.e-2 )   &
          write(debug%otun, *) " MMMMMLL"

!   output qni field:
        IF ( MINVAL(Cloud_state%qni_in+Cloud_state%SNi_out) .LT. -1.e-5) &
                                         write(Debug%otun, *) " MMMMMN2"

!   output delta T field:
        IF ( MAXVAL(ST_out) .GT. 7. ) write(Debug%otun, *) " MMMMMT1 "
        IF ( MINVAL(ST_out) .LT. - 7. ) write(Debug%otun, *) " MMMMMT2 "

!   output T field:
        IF ( MAXVAL(Input_mp%tin+ST_out) .GT. 330. )   &
                                          write(debug%otun, *) " MMMMMT3 "
        IF ( MINVAL(Input_mp%tin+ST_out) .LT. 170. )    &
                                          write(Debug%otun, *) " MMMMMT4 "

!   output 3D rain field:
        IF  ( MINVAL(Removal_mp%rain3d) .LT. 0. )   &
                                            write(Debug%otun, *) " MMMMMR1 "
!   output 3D snow field:
        IF  ( MINVAL(Removal_mp%snow3d) .LT. 0. )   &
                                           write(Debug%otun, *) " MMMMMS1 "

!   output sfc precip field:
        IF ( MINVAL(Precip_state%surfrain) .LT. 0. )   &
                                          write(Debug%otun, *) " MMMMMX1 "
!   output sfc snowfall field:
        IF ( MINVAL(Precip_state%surfsnow) .LT. 0. )   &
                                          write(Debug%otun, *) " MMMMMX2 "

      endif ! (debugo)

!-------------------------------------------------------------------------


end subroutine lscloud_debug 


!########################################################################

subroutine macrophysics_debug (Cloud_processes, Cloud_state)

!-------------------------------------------------------------------------
!    subroutine macrophysics_debug outputs requested debug info from the
!    model macrophysics modules.
!-------------------------------------------------------------------------

type(cloud_processes_type), intent(in) :: Cloud_processes
type(cloud_state_type),     intent(in) :: Cloud_state     


!--------------------------------------------------------------------------
!    output large-scale condensation, erosion and updated cloud area at
!    requested debug point.
!--------------------------------------------------------------------------
      write( Debug%otun,*) " dcond_ls ",    &
           Cloud_processes%dcond_ls(Debug%isamp,  Debug%jsamp, Debug%ksamp)
      write( Debug%otun,*) " D_eros ",   &
           Cloud_processes%D_eros(Debug%isamp,Debug%jsamp,Debug%ksamp)
      write( Debug%otun,*) " qa_upd ",  &
           Cloud_state%qa_upd(Debug%isamp,Debug%jsamp,Debug%ksamp)

!-----------------------------------------------------------------------

end subroutine macrophysics_debug 


!########################################################################

subroutine microphysics_debug1 (k, ST, label)

!-------------------------------------------------------------------------
!    subroutine microphysics_debug1 outputs temperature tendency from the
!    model microphysics module. 
!-------------------------------------------------------------------------

integer,                    intent(in) :: k
real, dimension(:,:,:),     intent(in) :: ST
character(len=*),           intent(in) :: label

!--------------------------------------------------------------------------
!    output temperature tendency at requested debug point.
!--------------------------------------------------------------------------
      if (k == Debug%ksamp) then
        write(Debug%otun,*) label,   ST(Debug%isamp,Debug%jsamp,Debug%ksamp)
      endif

!-----------------------------------------------------------------------

end subroutine microphysics_debug1 


!########################################################################

subroutine microphysics_debug2 (k, D_eros, D2_dt, tmp2, label)

!-------------------------------------------------------------------------
!    subroutine microphysics_debug2 outputs cloud erosion and terms
!    related to the temperature tendency due to liquid to ice conversions.
!-------------------------------------------------------------------------

integer,                    intent(in) :: k
real, dimension(:,:,:),     intent(in) :: D_eros, D2_dt, tmp2
character(len=*),           intent(in) :: label

!--------------------------------------------------------------------------
!    output requested fields at requested debug point.
!--------------------------------------------------------------------------
      if (k == Debug%ksamp) then
        write(Debug%otun,*) label,    &
                        D_eros(Debug%isamp, Debug%jsamp,Debug%ksamp),   &
                        D2_dt(Debug%isamp,Debug%jsamp,Debug%ksamp),&
                        tmp2(Debug%isamp,Debug%jsamp,Debug%ksamp)
      endif

!-----------------------------------------------------------------------

end subroutine microphysics_debug2 


!########################################################################

subroutine aerosol_cloud_debug1 (i, j, k, u_i, qs, qvsi, qvt, qin, &
                                 cf, cha, qa, chr) 

!-------------------------------------------------------------------------
!    subroutine aerosol_cloud_debug1 outputs relevant variables when the
!    relative humidity exceeds 200%.
!-------------------------------------------------------------------------

integer,      intent(in) :: i, j, k
real,         intent(in) :: u_i, qs, qvsi, qvt, qin, &
                            cf, cha, qa, chr

!--------------------------------------------------------------------------
!    output requested fields at requested debug point.
!--------------------------------------------------------------------------
      write(Debug%otun,*) " +++++++++++++++++++++++ "
      write(Debug%otun,*) " i,j,k ,u_i ", i,j,k,u_i
      write(Debug%otun,*) " qs, qvsi , qvt, qv ",  qs, qvsi, qvt, qin
      write(Debug%otun,*) " cf, ahuco,qa_upd,qrat ", cf,  cha, qa, chr
      write(Debug%otun,*) " qv(i,k) - cf * qs(i,k), (1-cf) " , &
                                                      qin - cf*qs,  (1-cf)
      write(Debug%otun,*) " +++++++++++++++++++++++ "

!-----------------------------------------------------------------------

end subroutine aerosol_cloud_debug1 


!########################################################################

subroutine aerosol_cloud_debug2 (i, j, k, rh)

!-------------------------------------------------------------------------
!    subroutine aerosol_cloud_debug1 outputs relevant variables when the
!    relatuive humidity exceeds 200%.
!-------------------------------------------------------------------------

integer,    intent(in) :: i,j,k
real,       intent(in) :: rh

!--------------------------------------------------------------------------
!    output requested fields at requested debug point.
!--------------------------------------------------------------------------
      write(Debug%otun,*) "MMMMMM RH", i,j,k
      write( Debug%otun,*) " rh_crit_min_1d ", rh   

!-----------------------------------------------------------------------

end subroutine aerosol_cloud_debug2 



!##########################################################################

subroutine write_debug_output (label, field, j, scalar_in)

!-------------------------------------------------------------------------
!    subroutine write_debug_output writes header info into the debug 
!    output file.
!-------------------------------------------------------------------------

character(len=*),       intent(in)           :: label
real, dimension(:,:,:), intent(in)           :: field
real,                   intent(in), optional :: scalar_in
integer,                intent(in), optional :: j

!----------------------------------------------------------------------
!   local variables:

      real :: scalar

!-----------------------------------------------------------------------
!    execute only if debug is active.
!-----------------------------------------------------------------------
      if (Debug%debugo) then

!-------------------------------------------------------------------------
!    define the value to use for scalar,  an optional factor which may be 
!    used to multiply the output field to put it in desired units before 
!    being written to the output file.
!-----------------------------------------------------------------------
        if (present(scalar_in) ) then
          scalar = scalar_in
        else
          scalar = 1.0
        endif

!-------------------------------------------------------------------------
!    if optional argument j is present, then only values for rows with 
!    j = jsamp, the supplied jindex for output, will be output. this 
!    argument is needed when the routine calling this debug routine is
!    being executed within a j loop, rather than in a single call over all
!    j indices of the processor / physics_window.
!-------------------------------------------------------------------------
        if (present(j) ) then
          if (j .eq. Debug%jsamp) then
            write (Debug%otun, *) label,    &
                        field(Debug%isamp, Debug%jsamp, Debug%ksamp)*scalar
          endif
        else
          write (Debug%otun, *) label,    &
                        field(Debug%isamp, Debug%jsamp, Debug%ksamp)*scalar
        endif
      endif
     
!-----------------------------------------------------------------------


end subroutine write_debug_output 



!########################################################################

subroutine record_micro_refusal

!-------------------------------------------------------------------------
!    record any cases where the microphysics was not executed because of 
!    negative water in the column.  this only an option with the "mg"
!    microphysics.
!-------------------------------------------------------------------------
      Debug%nrefuse = Debug%nrefuse + 1

!------------------------------------------------------------------------


end subroutine record_micro_refusal


!########################################################################

subroutine output_refusals 

!-------------------------------------------------------------------------
!    output to the stdout file a record of microphysics refusals. this 
!    only an option with "mg" microphysics.
!-------------------------------------------------------------------------

      if (Debug%debugo4 .and. Debug%nrefuse .gt. 0) then
        write(Debug%otun, '(a, i6, a, i6,a, / ,a, i6)')   &
             " in lscloud_debug_mod/output_refusals: pe = ", mpp_pe(),   &
                 "  timestep = ", Debug%ncall, " :",   &
                   "# of columns with negative water and thus skipping  &
                                            &microphysics: ", Debug%nrefuse
      end if

!-------------------------------------------------------------------------

end subroutine output_refusals


!########################################################################



                    end module lscloud_debug_mod
