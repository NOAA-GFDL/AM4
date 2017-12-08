module atmos_regional_tracer_driver_mod
!
! <CONTACT EMAIL="Larry.Horowitz@noaa.gov">
!   Larry W. Horowitz
! </CONTACT>

! <OVERVIEW>
!     This code calculates tendencies for regional chemical tracers
! </OVERVIEW>

! <DESCRIPTION>
!
! This code calculates emissions and chemical production/loss of regional
! chemical tracers tracers.
!
! Add more chemical tracers in extra regions like: central Africa, Russian Asia
! and southeastern Asia ---- by yyf, Jan 27, 2009
!
! </DESCRIPTION>


!-----------------------------------------------------------------------

use                    fms_mod, only : file_exist,   &
                                       field_exist, &
                                       write_version_number, &
                                       mpp_pe,  &
                                       mpp_root_pe, &
                                       open_namelist_file, &
                                       close_file,   &
                                       stdlog, &
                                       check_nml_error, &
                                       error_mesg, &
                                       FATAL, &
                                       WARNING
use           time_manager_mod, only : time_type
use           diag_manager_mod, only : send_data,            &
                                       register_diag_field
use         tracer_manager_mod, only : get_tracer_index,     &
                                       query_method
use          field_manager_mod, only : MODEL_ATMOS,          &
                                       parse
use              constants_mod, only : WTMAIR, AVOGNO
use           interpolator_mod, only : interpolate_type,     &
                                       interpolator_init,    &
                                       interpolator_end,     &
                                       interpolator,         &
                                       query_interpolator,   &
                                       CONSTANT,             &
                                       INTERP_WEIGHTED_P  

implicit none

private

!-----------------------------------------------------------------------
!     ... interfaces
!-----------------------------------------------------------------------
public  regional_tracer_driver, regional_tracer_driver_init


!-----------------------------------------------------------------------
!     ... namelist
!-----------------------------------------------------------------------
!  set a maximum number of tracers to assign the array
!  and update the actual number (ntracers) from the namelist
!  ntracers is read in regional_tracer_driver_init and 
!  must be publicly accessiable in the entire module
!-----------------------------------------------------------------------
integer            :: ntracers
integer, parameter :: max_ntracers = 50
real, dimension(max_ntracers) :: loss_freq = 4.63e-7
real                          :: loss_freq_all = 0.
real                          :: co_yield_from_avoc = 0.7, &
                                 co_yield_from_bvoc = 0.4, &
                                 co_yield_from_ch4  = 0.86,&
                                 tch4_vmr = 0.
character(len=16), dimension(max_ntracers) :: tracer_names = " "
 
namelist /regional_tracer_driver_nml/     &
                               ntracers,  &
                               loss_freq, &
                               loss_freq_all, &
                               tracer_names, &
                               co_yield_from_avoc, &
                               co_yield_from_bvoc, &
                               co_yield_from_ch4, &
                               tch4_vmr

!-----------------------------------------------------------------------
!     ...  declare type that will store the field infomation for the 
!          emission file
!-----------------------------------------------------------------------

type,public :: field_init_type
   character(len=64), pointer :: field_names(:)
end type field_init_type



character(len=7), parameter :: module_name = 'tracers'
real, parameter :: g_to_kg    = 1.e-3,    & !conversion factor (kg/g)
                   m2_to_cm2  = 1.e4        !conversion factor (cm2/m2)
real, parameter :: emis_cons = WTMAIR * g_to_kg * m2_to_cm2 / AVOGNO
logical, dimension(max_ntracers) :: Lemis = .false.
type(interpolate_type),dimension(max_ntracers) :: inter_emis
type(field_init_type),dimension(max_ntracers)  :: init_type
     
integer, dimension(max_ntracers) :: tracer_indices = 0

logical :: module_is_initialized=.false.

!-----------------------------------------------------------------------
!     ... identification numbers for diagnostic fields
!-----------------------------------------------------------------------
integer, dimension(max_ntracers) :: id_prod, id_loss, id_chem_tend, id_emiss

!-----------------------------------------------------------------------
!     ... tracer numbers
!-----------------------------------------------------------------------
integer :: id_tch4=0, id_avoc=0, id_bvoc=0, &
           id_cofromch4=0, id_cofromavoc=0, id_cofrombvoc=0

!---- version number ---------------------------------------------------
character(len=128), parameter :: version     = ''
character(len=128), parameter :: tagname     = ''
!-----------------------------------------------------------------------

contains


!#######################################################################

! <SUBROUTINE NAME="regional_tracer_driver">
!   <OVERVIEW>
!     Regional tracer driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calculates the sources and sinks of regional tracers
!     for "CO" and "SA"
!   </DESCRIPTION>
!   <TEMPLATE>
!     call regional_tracer_driver( lon, lat, pwt, r, chem_dt, &
!                                  Time, phalf, is, js, kbot)
!   </TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     Pressure weighting (air mass) for each layer (kg/m2)
!   </IN>
!   <IN NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (tropchem tracers in VMR)
!   </IN>
!   <IN NAME="Time, Time_next" TYPE="time_type">
!     Model time
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <IN NAME="is, js" TYPE="integer">
!     Local domain start indices
!   </IN>
!   <OUT NAME="chem_dt" TYPE="real" DIM="(:,:,:,:)">
!     Tracer tendencies from tropospheric chemistry (VMR/s)
!   </OUT>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

subroutine regional_tracer_driver( lon, lat, pwt, r, chem_dt, &
                                   Time, phalf, is, js, kbot)

!-----------------------------------------------------------------------
   real, intent(in),    dimension(:,:)            :: lon, lat
   real, intent(in),    dimension(:,:,:)          :: pwt
   real, intent(in),    dimension(:,:,:,:)        :: r
   real, intent(out),   dimension(:,:,:,:)        :: chem_dt
   type(time_type), intent(in)                    :: Time
   integer, intent(in)                            :: is,js
   real, intent(in),    dimension(:,:,:)          :: phalf
   integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(r,1),size(r,2)) :: emis
   real, dimension(size(r,1),size(r,2),ntracers) :: emis_save
   integer :: i,j,n,kb,id,jd,kd,trind
   integer :: nt, ntp
   logical :: used
   real, dimension(size(r,1),size(r,2),size(r,3),ntracers) :: emis_source
   character(len=64) :: name
   real, dimension(size(r,1),size(r,2),size(r,3),ntracers) :: prod, loss

   integer :: omp_get_num_threads
!-----------------------------------------------------------------------

! Not yet tested for OpenMP
!$ if(omp_get_num_threads() > 1) &
!$    call error_mesg ('Regional_tracer_driver','Code is not tested for OpenMP and omp_num_threads > 1', FATAL)

!<ERROR MSG="regional_tracer_driver_init must be called first." STATUS="FATAL">
!   Regional_tracer_driver_init needs to be called before regional_tracer_driver.
!</ERROR>
   if (.not. module_is_initialized)  &
      call error_mesg ('Regional_tracer_driver','regional_tracer_driver_init must be called first.', FATAL)

   id=size(r,1); jd=size(r,2); kd=size(r,3)
   nt = size(r,4); ntp = size(chem_dt,4)
 
   emis_source(:,:,:,:) = 0.0
   chem_dt(:,:,:,:) = 0.0

   do n = 1, ntracers
      trind = tracer_indices(n)

!-----------------------------------------------------------------------
!     ... read in the surface emissions, using interpolator
!-----------------------------------------------------------------------
      if (Lemis(n)) then
         call read_2D_emis_data( inter_emis(n), emis, Time, &
                                 init_type(n)%field_names, &
                                 is, js )
         emis_save(:,:,n) = emis(:,:)

         if (present(kbot)) then
            do j=1,jd
               do i=1,id
                  kb=kbot(i,j)
                  emis_source(i,j,kb,n) = emis(i,j)/pwt(i,j,kb) * emis_cons
                  
               end do
            end do
         else
            emis_source(:,:,kd,n) = emis(:,:)/pwt(:,:,kd) * emis_cons
         end if

      end if

!-----------------------------------------------------------------------
!     ... calculate chemical losses
!-----------------------------------------------------------------------
      if (trind > 0) then
         if( n==id_tch4 .and. tch4_vmr>0.) then
            loss(:,:,:,n) = loss_freq(n) * tch4_vmr
         else
            loss(:,:,:,n) = loss_freq(n) * r(:,:,:,trind)
         end if
      else
         loss(:,:,:,n) = 0.
      end if

   end do

!-----------------------------------------------------------------------
!     ... calculate chemical production
!-----------------------------------------------------------------------
   prod(:,:,:,:) = 0.
   if (id_tch4>0 .and. id_cofromch4>0) then
      prod(:,:,:,id_cofromch4)  = co_yield_from_ch4  * loss(:,:,:,id_tch4)
   end if
   if (id_avoc>0 .and. id_cofromavoc>0) then
      prod(:,:,:,id_cofromavoc) = co_yield_from_avoc * loss(:,:,:,id_avoc)
   end if
   if (id_bvoc>0 .and. id_cofrombvoc>0) then
      prod(:,:,:,id_cofrombvoc) = co_yield_from_bvoc * loss(:,:,:,id_bvoc)
   end if

   do n = 1, ntracers
      trind = tracer_indices(n)

!-----------------------------------------------------------------------
!     ... compute tendency
!-----------------------------------------------------------------------
      if (trind>0 .and. trind<=ntp) then
         chem_dt(:,:,:,trind) = emis_source(:,:,:,n) + prod(:,:,:,n) - loss(:,:,:,n)
      end if

!-----------------------------------------------------------------------
!     ... output diagnostics
!-----------------------------------------------------------------------
      if(id_emiss(n) > 0) then
         used = send_data(id_emiss(n), emis_save(:,:,n), Time, is_in=is, js_in=js)
      end if

      if(id_loss(n)>0) then
         used = send_data(id_loss(n),loss(:,:,:,n),Time,is_in=is,js_in=js)
      end if
      
      if(id_prod(n)>0) then
         used = send_data(id_prod(n),prod(:,:,:,n),Time,is_in=is,js_in=js)
      end if

      if(id_chem_tend(n)>0 .and. trind>0 .and. trind<=ntp) then
         used = send_data( id_chem_tend(n), chem_dt(:,:,:,trind),&
                     Time, is_in=is,js_in=js)
      end if

   end do     
   
!-----------------------------------------------------------------------
    
end subroutine regional_tracer_driver
!</SUBROUTINE>

!#######################################################################

! <FUNCTION NAME="regional_tracer_driver_init">
!   <OVERVIEW>
!     Initializes the regional tracer driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the regional tracer module.
!     It is called from atmos_tracer_driver_init.
!     Data sets are read in for surface emissions.
!   </DESCRIPTION>
!   <TEMPLATE>
!       call regional_tracer_driver_init (lonb_mod, latb_mod, axes, Time, mask)
!   </TEMPLATE>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask that designates which grid points
!      are above (1) or below (0) the ground
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="lonb_mod" TYPE="real" DIM="(:,:)">
!     The longitude corners for the local domain.
!   </IN>
!   <IN NAME="latb_mod" TYPE="real" DIM="(:,:)">
!     The latitude corners for the local domain.
!   </IN>

subroutine regional_tracer_driver_init( lonb_mod, latb_mod, axes, Time, mask )

!-----------------------------------------------------------------------
!
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
   type(time_type), intent(in) :: Time
   integer        , intent(in) :: axes(4)
   real, intent(in), dimension(:,:) :: lonb_mod
   real, intent(in), dimension(:,:) :: latb_mod
   real, intent(in),    dimension(:,:,:), optional :: mask

   integer :: n, trind
   character(len=64) :: trname
   character(len=64) :: diag_name
   integer :: ierr, io
   integer :: unit

   integer :: omp_get_num_threads


! Not yet tested for OpenMP
!$ if(omp_get_num_threads() > 1) &
!$    call error_mesg ('Regional_tracer_driver_init','Code is not tested for OpenMP and omp_num_threads > 1', FATAL)
!
!-----------------------------------------------------------------------
!
   if (module_is_initialized) return

!-----------------------------------------------------------------------
!     ... write version number
!-----------------------------------------------------------------------
   call write_version_number(version, tagname)
    
!-----------------------------------------------------------------------
!     ... read namelist
!-----------------------------------------------------------------------
   if(file_exist('input.nml')) then
      unit = open_namelist_file('input.nml')
      ierr=1; do while (ierr /= 0)
      read(unit, nml = regional_tracer_driver_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'regional_tracer_driver_nml')
      end do
10    call close_file(unit)
   end if
  
   if(mpp_pe() == mpp_root_pe()) then       
      write(stdlog(), nml=regional_tracer_driver_nml)
   end if
   
   
   if (loss_freq_all > 0.) then
      loss_freq(:) = loss_freq_all
   end if

!-----------------------------------------------------------------------
!     ... set initial value of indices
!-----------------------------------------------------------------------
   tracer_indices(:) = 0
   do n=1,ntracers
      trind = get_tracer_index(MODEL_ATMOS, tracer_names(n))
      if (trind >0) then
         tracer_indices(n) = trind
         write(*,30) tracer_names(n),tracer_indices(n)
      else
!<ERROR MSG="Tropospheric chemistry tracer not found in field table" STATUS="WARNING">
!   A tropospheric chemistry tracer was not included in the field table
!</ERROR>
         call error_mesg ('tropchem_driver_init', trim(tracer_names(n)) // ' is not found', WARNING)
      end if
   end do
30 format (A,' was initialized as tracer number ',i3)

!-----------------------------------------------------------------------
!     ... Keep track of tracer indices
!-----------------------------------------------------------------------
   do n=1,ntracers
      trname = tracer_names(n)
      select case (trname)
         case('tch4')
            id_tch4 = n
         case('avoc')
            id_avoc = n
         case('bvoc')
            id_bvoc = n
         case('cofromch4')
            id_cofromch4 = n
         case('cofromavoc')
            id_cofromavoc = n
         case('cofrombvoc')
            id_cofrombvoc = n
      end select
   end do

!-----------------------------------------------------------------------
!     ... Setup emissions input/interpolation
!-----------------------------------------------------------------------
   do n = 1,ntracers
      trind = tracer_indices(n)
      if (trind > 0) then
         call init_2D_emis_data( inter_emis(n), MODEL_ATMOS, tracer_indices(n),&
                                 lonb_mod, latb_mod, init_type(n), Lemis(n))
      end if
   end do        

!-----------------------------------------------------------------------
!     ... Register diagnostic fields for species tendencies
!-----------------------------------------------------------------------
   do n=1,ntracers
      trname = tracer_names(n)
      diag_name = trim(trname)//'_chem_dt'
      id_chem_tend(n) = register_diag_field( module_name, diag_name, axes(1:3), &
                                             Time, diag_name,'VMR/s' )
      diag_name = trim(trname)//'_loss'
      id_loss(n) = register_diag_field( module_name, diag_name, axes(1:3), &
                                        Time, diag_name,'VMR/s')
      diag_name = trim(trname)//'_prod'
      id_prod(n) = register_diag_field( module_name, diag_name, axes(1:3), &
                                        Time, diag_name,'VMR/s')
      if( Lemis(n) ) then
         diag_name = trim(trname)//'_emis'
         id_emiss(n) = register_diag_field( module_name, diag_name, axes(1:2), &
                                            Time, diag_name, 'molec/cm2/s')
      else
         id_emiss(n) = 0
      end if
   end do

   module_is_initialized = .true.
      
      
!-----------------------------------------------------------------------
      
end subroutine regional_tracer_driver_init
!</FUNCTION>
 
      
subroutine regional_tracer_driver_end

!-----------------------------------------------------------------------
!     ... initialize mpp clock id
!-----------------------------------------------------------------------
      
   module_is_initialized = .false.
      
      
!-----------------------------------------------------------------------
      
end subroutine regional_tracer_driver_end

!#######################################################################

! <SUBROUTINE NAME="read_2D_emis_data">
!   <OVERVIEW>
!     Read emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Reads tracer surface emissions from a NetCDF file
!   </DESCRIPTION>
!   <TEMPLATE>
!     call read_2D_emis_data( emis_type, emis, Time, field_names, is, js ) 
!   </TEMPLATE>

subroutine read_2D_emis_data( emis_type, emis, Time, field_names, is, js )
    
   type(interpolate_type),intent(inout) :: emis_type
   real, dimension(:,:),intent(out) :: emis
   type(time_type),intent(in) :: Time
   character(len=*),dimension(:), intent(in) :: field_names
   integer, intent(in) :: is, js


   integer :: k
   logical :: used
   real, dimension(size(emis,1),size(emis,2)) :: temp_data

   emis(:,:) = 0.
   temp_data(:,:) = 0.
   do k = 1,size(field_names)
      call interpolator(emis_type,Time,temp_data,field_names(k),is,js)
      emis(:,:) = emis(:,:) + temp_data(:,:)
   end do

end subroutine read_2D_emis_data
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="init_2D_emis_data">
!   <OVERVIEW>
!     Open emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Opens NetCDF file of tracer surface emissions for reading,
!     and set up interpolation to model grid/time
!   </DESCRIPTION>
!   <TEMPLATE>
!     call init_2D_emis_data( emis_type, model, index, &
!                             lonb_mod, latb_mod, field_type, flag )
!   </TEMPLATE>

subroutine init_2D_emis_data( emis_type, model, index, &
                              lonb_mod, latb_mod, field_type, flag )
    
   type(interpolate_type), intent(inout) :: emis_type
   integer,                intent(in)    :: model, index
   real, dimension(:,:),   intent(in)    :: lonb_mod, latb_mod
   type(field_init_type),  intent(out)   :: field_type
   logical,                intent(out)   :: flag
    
   character(len=64) :: name, control
   integer :: nfields
   integer :: flag_file
   character(len=64) :: emis_name, emis_file

   flag = .false.
   control = ''
   if( query_method('emissions',model,index,name,control) ) then
      if( trim(name) == 'file' ) then
         flag_file = parse(control, 'file', emis_file)
         if(flag_file > 0) then
            flag = .true.
            call interpolator_init( emis_type, trim(emis_file), &
                                    lonb_mod, latb_mod,  &
                                    data_out_of_bounds=(/CONSTANT/), &
                                    vert_interp=(/INTERP_WEIGHTED_P/) )
            call query_interpolator(emis_type,nfields=nfields)
            allocate(field_type%field_names(nfields))
            call query_interpolator(emis_type,field_names=field_type%field_names)
         end if
      end if
   end if
end subroutine init_2D_emis_data
!</SUBROUTINE>


!############################################################################
end module atmos_regional_tracer_driver_mod
