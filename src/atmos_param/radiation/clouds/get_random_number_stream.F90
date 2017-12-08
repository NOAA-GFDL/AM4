
module random_number_streams_mod

!----------------------------------------------------------------------

use mpp_mod,            only:  input_nml_file
use fms_mod,            only:  open_namelist_file, mpp_pe, &
                               mpp_root_pe, stdlog, fms_init, &
                               write_version_number, file_exist, &
                               check_nml_error, close_file, &
                               error_mesg, FATAL, NOTE
use time_manager_mod,   only:  time_type
use constants_mod,      only:  RADIAN

use random_numbers_mod, only:  randomNumberStream,   &
                               initializeRandomNumberStream, &
                               constructSeed

use cloudrad_types_mod, only:  cloudrad_control_type

!--------------------------------------------------------------------

implicit none 
private

!--------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'
!---------------------------------------------------------------------
!-------  interfaces --------

public  random_number_streams_init, &
        random_number_streams_end, &
        get_random_number_streams

!---------------------------------------------------------------------
!-------- namelist  ---------

logical :: force_use_of_temp_for_seed = .true.
                                        ! if true, when using stochastic 
                                        ! clouds, force the seed to use 
                                        ! top-model-level temps as input to
                                        ! random number generator
logical :: do_legacy_seed_generation = .false.
                                        ! setting this variable to .true.
                                        ! (not recommended except to
                                        ! reproduce previous results) will
                                        ! activate the seed generation
                                        ! scheme used previously (through
                                        ! the siena code). this scheme may
                                        ! exhibit flaws at hi-res or when
                                        ! time is held fixed in the 
                                        ! radiation calculation.

namelist /random_number_streams_nml/ do_legacy_seed_generation, &
                                     force_use_of_temp_for_seed

!---------------------------------------------------------------------
!----  private data -------

logical :: module_is_initialized = .false.

!----------------------------------------------------------------------
!   variables needed for legacy random number seed:
!----------------------------------------------------------------------
real, dimension(:,:), allocatable  :: lats, lons ! lat and lon of columns
                                                 ! in this processor's
                                                 ! domain [ degrees ]
logical :: use_temp_for_seed


CONTAINS

!######################################################################

subroutine random_number_streams_init ( lonb, latb, Cldrad_control )
real, dimension(:,:),        intent(in)    ::  lonb, latb
type(cloudrad_control_type), intent(inout) ::  Cldrad_control

!----------------------------------------------------------------------
!   local variables:

      integer  ::   unit, ierr, io
      integer  ::   id, jd, i, j, ii, jj

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init

!---------------------------------------------------------------------
!    read namelist.
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=random_number_streams_nml, iostat=io)
      ierr = check_nml_error(io,"random_number_streams_nml")
#else
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read (unit, nml=random_number_streams_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'random_number_streams_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
           write (stdlog(), nml=random_number_streams_nml)

!---------------------------------------------------------------------
!     determine the source of the random number seed generator to be used
!     for the stochastic cloud generation. the legacy scheme may fail to
!     provide spacially unique seeds at hi-res (above c48) or if the time
!     provided the radiation package does not monotonically advance (as in
!     some specialized sensitivity / assessment studies). 
!---------------------------------------------------------------------
      if (do_legacy_seed_generation) then
!---------------------------------------------------------------------
!     if it is desired to force the use of the temperature-based
!     random number seed (as is used when time is not always advancing
!     as seen by the radiation package, or when model resolution is
!     less than 1 degree), set the logical control variable in 
!     Cldrad_control to so indicate. 
!---------------------------------------------------------------------
        if ( force_use_of_temp_for_seed) then
          use_temp_for_seed = .true.
          call error_mesg ('cloud_spec_init', &
                 'Will use temp as basis for stochastic cloud seed; &
                    eed is set true', NOTE)
        else
          use_temp_for_seed = .false.
          call error_mesg ('cloud_spec_init', &
               ' If model resolution is above c48, it is &
               &HIGHLY RECOMMENDED that you set cloud_spec_nml variable &
               &force_use_of_temp_for_seed to true to assure &
               &reproducibility across pe count and domain layout', NOTE)
          call error_mesg ('cloud_spec_init', &
               'No action is needed at or below c48 resolution.', NOTE)
        endif

!---------------------------------------------------------------------
!     if the latitude and longitude of adjacent points on a pe have the 
!     same integral values (NINT), set the logical control variable in 
!     Cldrad_control to use the model temperature at the top level as the 
!     random number seed to provide spacial uniqueness at the points on the
!     processor. 
!     Note that for model resolutions of ~ 1 degree, some pes may use
!     lat and lon, while others use temperature, and that this may change
!     as the domain decomposition or npes used for the problem are changed.
!     Therefore it is  HIGHLY RECOMMENDED that for  resolutions above c48  
!     that nml variable force_use_of_temp_for_seed be set to .true.; 
!     it may remain the default value of .false. for lower resolution runs
!     or to preserve legacy results, or if reproducibility over npes or 
!     layout is not essential. 
!---------------------------------------------------------------------
        if (.not. use_temp_for_seed) then
          id = size(lonb,1) - 1
          jd = size(latb,2) - 1
!--------------------------------------------------------------------
!    allocate and define arrays holding the processor's latitudes
!    and longitudes
!--------------------------------------------------------------------
          allocate (lats(size(latb,1),size(latb,2)))
          allocate (lons(size(lonb,1), size(lonb,2)))
          lats(:,:) = latb(:,:)*RADIAN
          lons(:,:) = lonb(:,:)*RADIAN


  jLoop:  do j=1,jd
            do i=1,id
              do jj=j+1,jd+1
                do ii=i+1,id+1
                  if (NINT(lats(ii,jj)) == NINT(lats(i,j))) then
                    if (NINT(lons(ii,jj)) == NINT(lons(i,j))) then    
                      Cldrad_control%use_temp_for_seed = .true.
                      call error_mesg ('cloud_spec_init', &
                           'Found grid point within 1 degree of  &
                             &another',NOTE)
                      call error_mesg ('cloud_spec_init', &
                            'if reproducibility across npes and layout is &
                           &desired, you must set cloud_spec_nml variable &
                           &force_use_of_temp_for_seed to true., and &
                           &restart the model.', NOTE)
                      exit jLoop
                    endif
                  endif
                enddo
              enddo
            enddo
          enddo jLoop
        endif

!-------------------------------------------------------------------------
!    set seed generation source to be the temperature field (default)
!-------------------------------------------------------------------------
      else
        use_temp_for_seed = .true.
      endif

!---------------------------------------------------------------------
!    mark the module initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------

end subroutine random_number_streams_init
  
!######################################################################

subroutine get_random_number_streams ( is, js, Rad_time, temp, streams, perm )
integer,                  intent(in)                  :: is, js
type(time_type),          intent(in)                  :: Rad_time
real,                     intent(in),  dimension(:,:) :: temp
type(randomNumberStream), intent(out), dimension(:,:) :: streams
integer,                  intent(in),  optional       :: perm

integer :: i, j
real    :: seedwts(8) = (/3000.,1000.,300.,100.,30.,10.,3.,1./)

      if (use_temp_for_seed) then
        do j = 1, size(streams,2)
          do i = 1, size(streams,1)
            streams(i,j) = initializeRandomNumberStream(  &
                               ishftc(nint(temp(i,j)*seedwts),perm))

          enddo
        enddo
      else
        do j = 1, size(streams,2)
          do i = 1, size(streams,1)
            streams(i,j) = initializeRandomNumberStream(  &
                            constructSeed(nint(lons(is + i - 1, js + j - 1)), &
                                          nint(lats(is + i - 1, js + j - 1)), &
                                                    Rad_time, perm=perm))
          enddo
        enddo
      endif

!----------------------------------------------------------------------

end subroutine get_random_number_streams

!######################################################################

subroutine random_number_streams_end

!----------------------------------------------------------------------
!    skip if the module is not initialized.
!----------------------------------------------------------------------
      if (.not.module_is_initialized) return

!----------------------------------------------------------------------
!    deallocate module arrays
!----------------------------------------------------------------------
      if (allocated(lats)) deallocate(lats)
      if (allocated(lons)) deallocate(lons)

!----------------------------------------------------------------------
!    mark the module as no longer initialized.
!----------------------------------------------------------------------
      module_is_initialized = .false.

!----------------------------------------------------------------------

end subroutine random_number_streams_end

!######################################################################

end module random_number_streams_mod

