module solar_data_driver_mod

!   shared modules:

use mpp_mod,          only: input_nml_file
use fms_mod,          only: fms_init, mpp_pe, mpp_root_pe, &
                            open_namelist_file, stdlog, &
                            file_exist, FATAL, NOTE, &
                            error_mesg, close_file, &
                            write_version_number, check_nml_error
use time_manager_mod, only: time_type, get_date, time_manager_init

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!------------ version number for this module --------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

!---------------------------------------------------------------------
!------ interfaces -----

public solar_data_driver_init, &
       solar_data_driver_time_vary, &
       solar_data_driver_end

!---------------------------------------------------------------------
!------- namelist ---------

real    :: solar_scale_factor = -1.0  ! factor to multiply incoming solar 
                                      ! spectral irradiances. default is to
                                      ! perform no computation. used to 
                                      ! change "solar constant"

namelist /solar_data_driver_nml/ solar_scale_factor

!--------------------------------------------------------------------
! miscellaneous variables and indices

! total flux (solar constant)
real, dimension(:,:), allocatable :: solflxtot_lean
real                              :: solflxtot_lean_ann_1882, &
                                     solflxtot_lean_ann_2000
! flux by band
real, dimension(:,:,:), allocatable :: solflxband_lean
real, dimension(:),     allocatable :: solflxband_lean_ann_1882, &
                                       solflxband_lean_ann_2000

integer ::  first_yr_lean, last_yr_lean,   &
            nvalues_per_year_lean, numbands_lean
integer ::  years_of_data_lean = 0
logical ::  module_is_initialized = .false.


CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!######################################################################

subroutine solar_data_driver_init (nbands, ierror)

integer, intent(in)  :: nbands
integer, intent(out) :: ierror

!---------------------------------------------------------------------
!   local variables

integer :: yr, month, nband, nyr, nv
integer ::   unit, io, ierr, logunit

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!--------------------------------------------------------------------
      call fms_init
      call time_manager_init
!     call esfsw_parameters_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=solar_data_driver_nml, iostat=io)
      ierr = check_nml_error(io,'solar_data_driver_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=solar_data_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'solar_data_driver_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!--------------------------------------------------------------------
!     code to handle time-varying solar input
!--------------------------------------------------------------------
        if (file_exist('INPUT/lean_solar_spectral_data.dat')) then
          unit = open_namelist_file   &
                                 ('INPUT/lean_solar_spectral_data.dat')
          read (unit, FMT = '(4i8)') first_yr_lean, last_yr_lean,  &
                                   nvalues_per_year_lean, numbands_lean
          if (numbands_lean /= nbands) then
            call error_mesg ('radiation_driver_mod', &
            ' number of sw parameterization bands in solar_spectral &
            &data file differs from that defined in esfsw_parameters',&
                                                           FATAL)
          endif

          years_of_data_lean = last_yr_lean - first_yr_lean + 1
          allocate (solflxtot_lean (years_of_data_lean, nvalues_per_year_lean))
          allocate (solflxband_lean(years_of_data_lean, nvalues_per_year_lean, numbands_lean))
          allocate (solflxband_lean_ann_1882(numbands_lean))
          allocate (solflxband_lean_ann_2000(numbands_lean))

          read (unit, FMT = '(2i6,f17.4)') yr, month, solflxtot_lean_ann_1882
          read (unit, FMT = '(6e12.5 )')  (solflxband_lean_ann_1882(nband), nband =1,numbands_lean)
          if (solar_scale_factor >= 0.0) then
            solflxtot_lean_ann_1882 = solflxtot_lean_ann_1882 * solar_scale_factor
            do nband = 1,numbands_lean
               solflxband_lean_ann_1882(nband) =    &
                          solflxband_lean_ann_1882(nband) * solar_scale_factor
            enddo
          endif

          do nyr=1,years_of_data_lean
            do nv=1,nvalues_per_year_lean
              read (unit, FMT = '(2i6,f17.4)') yr, month, solflxtot_lean(nyr,nv)
              read (unit, FMT = '(6e12.5 )')  (solflxband_lean  &
                                (nyr,nv,nband), nband =1,numbands_lean)
            end do
          end do
          if (solar_scale_factor >= 0.0) then
            do nv=1,nvalues_per_year_lean
              do nyr=1,years_of_data_lean
                solflxtot_lean(nyr,nv) =    &
                         solflxtot_lean(nyr,nv) * solar_scale_factor
              enddo
            enddo
            do nband = 1,numbands_lean
              do nv=1,nvalues_per_year_lean
                do nyr=1,years_of_data_lean
                  solflxband_lean(nyr,nv,nband) =    &
                           solflxband_lean(nyr,nv,nband) * solar_scale_factor
                enddo
              enddo
            enddo
          endif

          read (unit, FMT = '(2i6,f17.4)') yr, month, solflxtot_lean_ann_2000
          read (unit, FMT = '(6e12.5 )')  (solflxband_lean_ann_2000(nband), nband =1,numbands_lean)
          call close_file (unit)
          if (solar_scale_factor >= 0.0) then
            solflxtot_lean_ann_2000 = solflxtot_lean_ann_2000 * solar_scale_factor
            do nband = 1,numbands_lean
              solflxband_lean_ann_2000(nband) =    &
                        solflxband_lean_ann_2000(nband) * solar_scale_factor
            enddo
          endif

          ierror = 0
        else
          ierror = 1
        endif

!----------------------------------------------------------------------

end subroutine solar_data_driver_init

!######################################################################

subroutine solar_data_driver_time_vary (Solar_time, solar_constant, solflxband)

type(time_type),    intent(in)  :: Solar_time
real,               intent(out) :: solar_constant
real, dimension(:), intent(out) :: solflxband

integer :: year, month, dum
!---------------------------------------------------------------------
!    define the solar_constant appropriate at Solar_time
!--------------------------------------------------------------------
!    define time to be used for solar input data.
!--------------------------------------------------------------------
      call get_date (Solar_time, year, month, dum, dum, dum, dum)

!--------------------------------------------------------------------
!    define input value based on year and month of Solar_time.
!--------------------------------------------------------------------
      call get_solar_data (year, month, solar_constant, solflxband)

!--------------------------------------------------------------------

end subroutine solar_data_driver_time_vary

!######################################################################

subroutine get_solar_data (year, month, solar_constant, solflxband)

integer,            intent(in)  :: year, month
real,               intent(out) :: solar_constant
real, dimension(:), intent(out) :: solflxband

integer :: nband
!--------------------------------------------------------------------
!    returns solar constant for year and month
!--------------------------------------------------------------------

!   first some error checks

      if (years_of_data_lean == 0) then
          call error_mesg ('solar_data_driver_mod', &
             'no lean data found, data file may not exist', FATAL)
      endif

      if (size(solflxband(:)) /= numbands_lean) then
        call error_mesg ('solar_data_driver_mod', &
              'bands present in solar constant time data differs from &
               &model parameterization band number', FATAL)
      endif

!--------------------------------------------------------------------

      if (year < first_yr_lean) then
         solar_constant = solflxtot_lean_ann_1882
         do nband=1,numbands_lean
            solflxband(nband) = solflxband_lean_ann_1882(nband)
         end do
      else if (year > last_yr_lean) then
         solar_constant = solflxtot_lean_ann_2000
         do nband=1,numbands_lean
            solflxband(nband) = solflxband_lean_ann_2000(nband)
         end do
      else
         solar_constant = solflxtot_lean(year-first_yr_lean+1, month)
         do nband=1,numbands_lean
            solflxband(nband) = solflxband_lean(year-first_yr_lean+1, month, nband)
         end do
      endif

!--------------------------------------------------------------------

end subroutine get_solar_data

!######################################################################

subroutine solar_data_driver_end

!--------------------------------------------------------------------

      deallocate (solflxtot_lean)
      deallocate (solflxband_lean)
      deallocate (solflxband_lean_ann_1882)
      deallocate (solflxband_lean_ann_2000)

      module_is_initialized = .false.

!--------------------------------------------------------------------
end subroutine solar_data_driver_end

!######################################################################

end module solar_data_driver_mod
