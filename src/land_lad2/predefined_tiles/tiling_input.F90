module predefined_tiles_mod

 use hdf5
 use,intrinsic :: iso_c_binding
 use constants_mod     , only : pi
 use land_data_mod, only: land_state_type
 use vegn_cohort_mod, only : vegn_cohort_type
 use land_tile_mod, only : first_elmt, insert, new_land_tile_predefined
 use land_tile_mod, only : land_tile_list_type,land_tile_type,&
                           land_tile_enum_type
 use tiling_input_types_mod, only : tile_parameters_type,lake_predefined_type
 use tiling_input_types_mod, only : glacier_predefined_type,soil_predefined_type
 use tiling_input_types_mod, only : metadata_predefined_type
 use time_manager_mod, only : time_type
 use time_interp_mod, only : time_interp

 implicit none

 interface get_parameter_data
   module procedure get_parameter_data_1d_integer
   module procedure get_parameter_data_1d_real
   module procedure get_parameter_data_2d_integer
   module procedure get_parameter_data_2d_real
 end interface

contains

subroutine open_database_predefined_tiles(h5id)

  integer :: status,h5id
  !Initialize the fortran library
  call h5open_f(status)

  !Open access to the model input database
  CALL h5fopen_f('INPUT/land_model_input_database.h5',H5F_ACC_RDONLY_F,h5id, status)
  !CALL h5fopen_f('INPUT/land_model_input_database.nc',H5F_ACC_RDONLY_F,h5id, status)

end subroutine open_database_predefined_tiles

subroutine close_database_predefined_tiles(h5id)

  integer :: status,h5id
  !Close access to the model input database
  call h5fclose_f(h5id,status)

  !Close the hdf5 library
  call h5close_f(status)

end subroutine close_database_predefined_tiles

subroutine determine_cell_id(is,js,h5id,cellid)
 
  integer,intent(in) :: is,js,h5id
  integer,intent(inout) :: cellid
  integer :: status,varid,dsid,grpid
  !integer :: dsid,cid
  real*8 :: h5tmp(1,1)
  integer(hsize_t) :: dims(2),maxdims(2)
  integer :: memrank
  integer(hsize_t) :: dimsm(2),count(2),offset(2)
  integer(hid_t) :: memid
 
 !NOTE: SHOULD BE DONE BY I/O CORE
 !Determine the id that corresponds to the lat/lon
 call h5gopen_f(h5id,"metadata",grpid,status)
 call h5dopen_f(grpid,"mapping",varid,status)
 !Get dataset's dataspace identifier.
 call h5dget_space_f(varid,dsid,status)
 ! Select hyperslab in the dataset
 offset = (/js-1,is-1/)
 count = (/1,1/)
 call h5sselect_hyperslab_f(dsid,H5S_SELECT_SET_F,offset,count,status)
 ! Create memory dataspace.
 memrank = 0
 dimsm = (/1,1/)
 call h5screate_simple_f(memrank,dimsm,memid,status)
 ! Select hyperslab in memory
 call h5sselect_hyperslab_f(memid,H5S_SELECT_SET_F,offset,count,status)
 ! Read data from hyperslab to memory
 dims = (/1,1/)
 call h5dread_f(varid,H5T_IEEE_F64LE,h5tmp,dims,status,memid,dsid)
 cellid = int(h5tmp(1,1))
 !print*,cellid

 !Close the memory space, dataset, and group
 call h5sclose_f(memid,status)
 call h5dclose_f(varid,status)
 call h5gclose_f(grpid,status)

end subroutine determine_cell_id

subroutine load_group_into_memory(tile,is,js,h5id,buf_ptr,buf_len,image_ptr)

  integer,intent(in) :: tile,is,js,h5id
  type(c_ptr),intent(inout) :: buf_ptr
  integer(size_t),intent(inout) :: buf_len
  character(kind=c_char),intent(inout),allocatable,dimension(:),target :: image_ptr
  integer :: status,grpid,cell_grpid,dstid
  !character(100) :: cellid_string
  character(100) :: tile_string,is_string,js_string,cellid_string
  integer(hid_t) :: fapl
  integer(size_t), parameter :: memory_increment = 1000000

 !Open access to the group in the database that contains all the group information
 call h5gopen_f(h5id,"grid_data",grpid,status)
 !Write the cell id to string
 !write(cellid_string,'(I10)') cellid
 !cellid_string = trim('g' // trim(adjustl(cellid_string)))
 write(tile_string,'(I10)') tile 
 write(is_string,'(I10)') is 
 write(js_string,'(I10)') js
 cellid_string = trim('tile:' // trim(adjustl(tile_string)) // ',is:' // &
                 trim(adjustl(js_string)) // ',js:' // trim(adjustl(is_string)))
 !The goal here is to load the desired group of the cellid into memory. This buffer 
 !will then be sent to the land model core. However, there is no direct way to do this
 !with a group instead we have to:
 !1.Create a new file in memory
 !2.Copy the desired group to this new file
 !3.Use the HDF5 api to then load this new file as a buffer
 !Ensure that we are always working in memory
 call h5pcreate_f(H5P_FILE_ACCESS_F,fapl,status)
 !Setting the third parameter to false ensures that we never write this file to disk
 call h5pset_fapl_core_f(fapl,memory_increment,.False.,status)
 !Although we create this file it is always in memory. It never gets written to disk
 call h5fcreate_f("buffer.hdf5",H5F_ACC_TRUNC_F,dstid,status,access_prp=fapl)
 !Close access to the property list
 call h5pclose_f(fapl,status)
 !Copy the group from the original database to the new file
 call h5ocopy_f(grpid,cellid_string,dstid,'data',status)
 !print*,status,cellid_string
 if (status .eq. -1)then
  print*,'This group does not exist in the database'
  stop
 endif
 !Flush the file
 call h5fflush_f(dstid,H5F_SCOPE_GLOBAL_F,status)
 !Determine the size of the desired group (which is now a file in memory...)
 buf_len = 0
 buf_ptr = C_NULL_PTR
 call h5fget_file_image_f(dstid,buf_ptr,int(0,size_t),status,buf_len)
 !Load the entire new file (i.e., desired group) into memory
 allocate(image_ptr(1:buf_len))
 buf_ptr = c_loc(image_ptr(1)(1:1))
 call h5fget_file_image_f(dstid,buf_ptr,buf_len,status)
 !Close the copied file (release memory)
 call h5fclose_f(dstid,status)
 !Close access to the grid data group
 call h5gclose_f(grpid,status)

end subroutine load_group_into_memory

subroutine open_image_file(buf_ptr,buf_len,image_ptr,dstid)

 integer,intent(inout) :: dstid
 type(c_ptr),intent(inout) :: buf_ptr
 integer(size_t),intent(inout) :: buf_len
 character(kind=c_char),intent(inout),allocatable,dimension(:),target :: image_ptr
 integer(hid_t) :: fapl
 integer :: status

 !This will all be don eon the land model side. The only purpose is to open access to
 !the buffer as if it is a file
 !Open a new fapl
 call h5pcreate_f(H5P_FILE_ACCESS_F,fapl,status)
 call h5pset_fapl_core_f(fapl,buf_len,.False.,status)
 !Assign the buffer to memory
 call h5pset_file_image_f(fapl,buf_ptr,buf_len,status)
 !Discard the buffer
 deallocate(image_ptr)
 buf_ptr = C_NULL_PTR
 !Open the file from memory
 dstid = 0
 call h5fopen_f("test_image",H5F_ACC_RDONLY_F,dstid,status,fapl)

end subroutine open_image_file

subroutine land_cover_cold_start_0d_predefined_tiles(tiles,lnd,l,h5id)
  
  type(land_tile_list_type),intent(inout) :: tiles
  type(land_state_type),intent(in) :: lnd
  integer,intent(in) :: l,h5id
  type(land_tile_type), pointer :: tile
  integer :: itile,tid,is,js
  !integer :: parent_id = 0
  integer :: status,varid,grpid,dimid,cell_grpid,cellid,dstid
  integer :: dsid,cid
  !real*8 :: h5tmp(1,1)
  !integer(hsize_t) :: dims(2),maxdims(2)
  !character(100) :: cellid_string
  real :: lat,lon,t0,t1
  real,allocatable,dimension(:) :: tmp
  type(tile_parameters_type) :: tile_parameters
  !integer :: memrank 
  !integer(hsize_t) :: dimsm(2),count(2),offset(2)
  !integer(hid_t) :: memid,fapl
  type(c_ptr) :: buf_ptr
  integer(size_t) :: buf_len
  character(kind=c_char),allocatable,dimension(:),target :: image_ptr

  !Determine the lat/lon of the grid cell (degrees)
  is = lnd%i_index(l)
  js = lnd%j_index(l)
  lon = 180.0*lnd%lon(l)/pi
  lat = 180.0*lnd%lat(l)/pi

  !Determine the cell id (I/O core)
  !call determine_cell_id(is,js,h5id,cellid)

  !Retrieve buffer and buffer length of desired group (I/O core)
  !call load_group_into_memory(cellid,h5id,buf_ptr,buf_len,image_ptr)
  call load_group_into_memory(lnd%face,is,js,h5id,buf_ptr,buf_len,image_ptr)

  !Use buffer and buffer length to open new image file (Land model core)
  call open_image_file(buf_ptr,buf_len,image_ptr,dstid)

  !Retrieve the group
  call h5gopen_f(dstid,'data',cid,status)

  !Metadata
  call retrieve_metadata(tile_parameters,cid)

  !Soil parameters
  call retrieve_soil_parameters(tile_parameters,cid)

  !Lake parameters
  call retrieve_lake_parameters(tile_parameters,cid)

  !Glacier parameters
  call retrieve_glacier_parameters(tile_parameters,cid)

  !Vegetation parameters

  !Create the tiles
  do itile = 1,tile_parameters%metadata%ntile
   !if (tile_parameters%metadata%frac(itile) .eq. 0.0)cycle
   tid = tile_parameters%metadata%tid(itile)
   select case (tile_parameters%metadata%ttype(itile))
    case(1)
     tile => new_land_tile_predefined(frac=tile_parameters%metadata%frac(itile),&
            glac=tid,glacier_predefined=tile_parameters%glacier,&
            itile=tid)
    case(2)
     tile => new_land_tile_predefined(frac=tile_parameters%metadata%frac(itile),&
            lake=tid,lake_predefined=tile_parameters%lake,&
            itile=tid)
    case(3)
     tile => new_land_tile_predefined(frac=tile_parameters%metadata%frac(itile),&
           soil=1,vegn=tile_parameters%soil%vegn(tid),&
           htag_j=tile_parameters%soil%hidx_j(tid),&
           htag_k=tile_parameters%soil%hidx_k(tid),&
           soil_predefined=tile_parameters%soil,itile=tid)
   end select
   call insert(tile,tiles)
  enddo

  !Close access to the grid cell's group
  call h5gclose_f(cid,status)
  call h5fclose_f(dstid,status)
  
end subroutine

subroutine retrieve_metadata(tile_parameters,cid)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer,intent(inout) :: cid
  type(metadata_predefined_type),pointer :: metadata
  integer :: dimid,grpid,ntile,nband,status,dsid,varid
  integer(hsize_t) :: dims(1),maxdims(1)
  allocate(tile_parameters%metadata)
  metadata => tile_parameters%metadata

  !Retrieve the group id
  call h5gopen_f(cid,"metadata",grpid,status)

  !Retrieve the number of tiles
  call h5dopen_f(grpid,"tile",varid,status)
  call h5dget_space_f(varid,dsid,status)
  call h5sget_simple_extent_dims_f(dsid,dims,maxdims,status)
  call h5dclose_f(varid,status)
  metadata%ntile = dims(1)

  !Retrieve the info
  call get_parameter_data(grpid,'frac',metadata%ntile,metadata%frac)
  call get_parameter_data(grpid,'tile',metadata%ntile,metadata%tile)
  call get_parameter_data(grpid,'type',metadata%ntile,metadata%ttype)
  call get_parameter_data(grpid,'tid',metadata%ntile,metadata%tid)

  !Clean up (This should go in the database creation)
  where ((metadata%frac .lt. 1.e-8) .and. (metadata%ttype .eq. 2))metadata%frac = 0.0
  metadata%frac = metadata%frac/sum(metadata%frac)

  !Close access to the group
  call h5gclose_f(grpid,status)

end subroutine retrieve_metadata

subroutine retrieve_glacier_parameters(tile_parameters,cid)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer,intent(inout) :: cid
  type(glacier_predefined_type),pointer :: glacier
  integer :: dimid,grpid,nglacier,nband,status,dsid,varid
  integer(hsize_t) :: dims(2),maxdims(2)
  allocate(tile_parameters%glacier)
  glacier => tile_parameters%glacier

  !Retrieve the group id
  call h5gopen_f(cid,"glacier",grpid,status)

  !If it does not exist then exit
  if (status .eq. -1)return

  !Retrieve the number of lake tiles and bands
  call h5dopen_f(grpid,"refl_min_dir",varid,status)
  call h5dget_space_f(varid,dsid,status)
  call h5sget_simple_extent_dims_f(dsid,dims,maxdims,status)
  nglacier = dims(1)
  nband = dims(2)
  glacier%nglacier = nglacier
  call h5sclose_f(dsid,status)
  call h5dclose_f(varid,status)

  !Retrieve the parameters
  call get_parameter_data(grpid,"w_sat",nglacier,glacier%w_sat)
  call get_parameter_data(grpid,"awc_lm2",nglacier,glacier%awc_lm2)
  call get_parameter_data(grpid,"k_sat_ref",nglacier,glacier%k_sat_ref)
  call get_parameter_data(grpid,"psi_sat_ref",nglacier,glacier%psi_sat_ref)
  call get_parameter_data(grpid,"chb",nglacier,glacier%chb)
  call get_parameter_data(grpid,"alpha",nglacier,glacier%alpha)
  call get_parameter_data(grpid,"heat_capacity_ref",nglacier,glacier%heat_capacity_ref)
  call get_parameter_data(grpid,"thermal_cond_ref",nglacier,glacier%thermal_cond_ref)
  call get_parameter_data(grpid,"emis_dry",nglacier,glacier%emis_dry)
  call get_parameter_data(grpid,"emis_sat",nglacier,glacier%emis_sat)
  call get_parameter_data(grpid,"z0_momentum",nglacier,glacier%z0_momentum)
  call get_parameter_data(grpid,"tfreeze",nglacier,glacier%tfreeze)
  call get_parameter_data(grpid,"refl_max_dir",nglacier,nband,glacier%refl_max_dir)
  call get_parameter_data(grpid,"refl_max_dif",nglacier,nband,glacier%refl_max_dif)
  call get_parameter_data(grpid,"refl_min_dir",nglacier,nband,glacier%refl_min_dir)
  call get_parameter_data(grpid,"refl_min_dif",nglacier,nband,glacier%refl_min_dif)

  !Close access to the group
  call h5gclose_f(grpid,status)

end subroutine retrieve_glacier_parameters

subroutine retrieve_lake_parameters(tile_parameters,cid)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer,intent(inout) :: cid
  type(lake_predefined_type),pointer :: lake
  integer :: dimid,grpid,nlake,nband,status,dsid,varid
  integer(hsize_t) :: dims(2),maxdims(2)
  allocate(tile_parameters%lake)
  lake => tile_parameters%lake

  !Retrieve the group id
  call h5gopen_f(cid,"lake",grpid,status)

  !If it does not exist then exit
  if (status .eq. -1)return

  !Retrieve the number of lake tiles and bands
  call h5dopen_f(grpid,"refl_dry_dir",varid,status)
  call h5dget_space_f(varid,dsid,status)
  call h5sget_simple_extent_dims_f(dsid,dims,maxdims,status)
  nlake = dims(1)
  nband = dims(2)
  lake%nlake = nlake
  call h5sclose_f(dsid,status)
  call h5dclose_f(varid,status)

  !Retrieve the parameters
  call get_parameter_data(grpid,"connected_to_next",nlake,lake%connected_to_next)
  call get_parameter_data(grpid,"whole_lake_area",nlake,lake%whole_area)
  call get_parameter_data(grpid,"lake_depth_sill",nlake,lake%depth_sill)
  call get_parameter_data(grpid,"lake_width_sill",nlake,lake%width_sill)
  call get_parameter_data(grpid,"lake_backwater",nlake,lake%backwater)
  call get_parameter_data(grpid,"lake_backwater_1",nlake,lake%backwater_1)
  call get_parameter_data(grpid,"awc_lm2",nlake,lake%awc_lm2)
  call get_parameter_data(grpid,"w_sat",nlake,lake%w_sat)
  call get_parameter_data(grpid,"chb",nlake,lake%chb)
  call get_parameter_data(grpid,"refl_dry_dif",nlake,nband,lake%refl_dry_dif)
  call get_parameter_data(grpid,"refl_dry_dir",nlake,nband,lake%refl_dry_dir)
  call get_parameter_data(grpid,"refl_sat_dif",nlake,nband,lake%refl_sat_dif)
  call get_parameter_data(grpid,"refl_sat_dir",nlake,nband,lake%refl_sat_dir)
  call get_parameter_data(grpid,"z0_momentum",nlake,lake%z0_momentum)
  call get_parameter_data(grpid,"psi_sat_ref",nlake,lake%psi_sat_ref)
  call get_parameter_data(grpid,"emis_dry",nlake,lake%emis_dry)
  call get_parameter_data(grpid,"z0_momentum_ice",nlake,lake%z0_momentum_ice)
  call get_parameter_data(grpid,"heat_capacity_ref",nlake,lake%heat_capacity_ref)
  call get_parameter_data(grpid,"alpha",nlake,lake%alpha)
  call get_parameter_data(grpid,"thermal_cond_ref",nlake,lake%thermal_cond_ref)
  call get_parameter_data(grpid,"emis_sat",nlake,lake%emis_sat)
  call get_parameter_data(grpid,"k_sat_ref",nlake,lake%k_sat_ref)

  !Close access to the group
  call h5gclose_f(grpid,status)

end subroutine retrieve_lake_parameters

subroutine retrieve_soil_parameters(tile_parameters,cid)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer,intent(inout) :: cid
  type(soil_predefined_type),pointer :: soil
  integer :: varid,grpid,nsoil,nband,status,dsid
  integer(hsize_t) :: dims(2),maxdims(2)
  allocate(tile_parameters%soil)
  soil => tile_parameters%soil

  !Retrieve the group id
  call h5gopen_f(cid,"soil",grpid,status)

  !If it does not exist then exit
  if (status .eq. -1)return

  !Retrieve the number of soil tiles and bands
  call h5dopen_f(grpid,"dat_refl_dry_dir",varid,status)
  call h5dget_space_f(varid,dsid,status)
  call h5sget_simple_extent_dims_f(dsid,dims,maxdims,status)
  nsoil = dims(1)
  nband = dims(2)
  soil%nsoil = nsoil
  call h5sclose_f(dsid,status)
  call h5dclose_f(varid,status)

  !Retrieve the parameters
  call get_parameter_data(grpid,"dat_w_sat",nsoil,soil%dat_w_sat)
  call get_parameter_data(grpid,"dat_awc_lm2",nsoil,soil%dat_awc_lm2)
  call get_parameter_data(grpid,"dat_k_sat_ref",nsoil,soil%dat_k_sat_ref)
  call get_parameter_data(grpid,"dat_psi_sat_ref",nsoil,soil%dat_psi_sat_ref)
  call get_parameter_data(grpid,"dat_chb",nsoil,soil%dat_chb)
  call get_parameter_data(grpid,"dat_heat_capacity_dry",nsoil,soil%dat_heat_capacity_dry)
  call get_parameter_data(grpid,"dat_thermal_cond_dry",nsoil,soil%dat_thermal_cond_dry)
  call get_parameter_data(grpid,"dat_thermal_cond_sat",nsoil,soil%dat_thermal_cond_sat)
  call get_parameter_data(grpid,"dat_thermal_cond_exp",nsoil,soil%dat_thermal_cond_exp)
  call get_parameter_data(grpid,"dat_thermal_cond_scale",nsoil,soil%dat_thermal_cond_scale)
  call get_parameter_data(grpid,"dat_thermal_cond_weight",nsoil,soil%dat_thermal_cond_weight)
  call get_parameter_data(grpid,"dat_refl_dry_dir",nsoil,nband,soil%dat_refl_dry_dir)
  call get_parameter_data(grpid,"dat_refl_dry_dif",nsoil,nband,soil%dat_refl_dry_dif)
  call get_parameter_data(grpid,"dat_refl_sat_dir",nsoil,nband,soil%dat_refl_sat_dir)
  call get_parameter_data(grpid,"dat_refl_sat_dif",nsoil,nband,soil%dat_refl_sat_dif)
  call get_parameter_data(grpid,"dat_emis_dry",nsoil,soil%dat_emis_dry)
  call get_parameter_data(grpid,"dat_emis_sat",nsoil,soil%dat_emis_sat)
  call get_parameter_data(grpid,"dat_z0_momentum",nsoil,soil%dat_z0_momentum)
  call get_parameter_data(grpid,"dat_tf_depr",nsoil,soil%dat_tf_depr)
  call get_parameter_data(grpid,"rsa_exp_global",nsoil,soil%rsa_exp_global)
  call get_parameter_data(grpid,"gw_res_time",nsoil,soil%gw_res_time)
  call get_parameter_data(grpid,"gw_hillslope_length",nsoil,soil%gw_hillslope_length)
  call get_parameter_data(grpid,"gw_scale_length",nsoil,soil%gw_scale_length)
  call get_parameter_data(grpid,"gw_hillslope_zeta_bar",nsoil,soil%gw_hillslope_zeta_bar)
  call get_parameter_data(grpid,"gw_hillslope_relief",nsoil,soil%gw_hillslope_relief)
  call get_parameter_data(grpid,"gw_scale_relief",nsoil,soil%gw_scale_relief)
  call get_parameter_data(grpid,"gw_soil_e_depth",nsoil,soil%gw_soil_e_depth)
  call get_parameter_data(grpid,"gw_scale_soil_depth",nsoil,soil%gw_scale_soil_depth)
  !call get_parameter_data(grpid,"gw_hillslope_a",nsoil,soil%gw_hillslope_a)
  !call get_parameter_data(grpid,"gw_hillslope_n",nsoil,soil%gw_hillslope_n)
  call get_parameter_data(grpid,"gw_perm",nsoil,soil%gw_perm)
  call get_parameter_data(grpid,"gw_scale_perm",nsoil,soil%gw_scale_perm)
  call get_parameter_data(grpid,"gw_res_time",nsoil,soil%gw_res_time)
  call get_parameter_data(grpid,"microtopo",nsoil,soil%microtopo)
  call get_parameter_data(grpid,"tile_hlsp_length",nsoil,soil%tile_hlsp_length)
  call get_parameter_data(grpid,"tile_hlsp_slope",nsoil,soil%tile_hlsp_slope)
  call get_parameter_data(grpid,"tile_hlsp_elev",nsoil,soil%tile_hlsp_elev)
  call get_parameter_data(grpid,"tile_hlsp_hpos",nsoil,soil%tile_hlsp_hpos)
  call get_parameter_data(grpid,"tile_hlsp_width",nsoil,soil%tile_hlsp_width)
  call get_parameter_data(grpid,"hidx_k",nsoil,soil%hidx_k)
  call get_parameter_data(grpid,"hidx_j",nsoil,soil%hidx_j)
  call get_parameter_data(grpid,"vegn",nsoil,soil%vegn)

  !Close access to the group
  call h5gclose_f(grpid,status)

end subroutine retrieve_soil_parameters

subroutine get_parameter_data_1d_integer(grpid,var,nx,tmp)

 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx
 integer,dimension(:),pointer :: tmp
 integer :: itile,varid,status
 integer(hsize_t) :: dims(1)
 allocate(tmp(nx))

 dims(1) = nx
 call h5dopen_f(grpid,var,varid,status)
 call h5dread_f(varid,H5T_STD_I32LE,tmp,dims,status)
 !print*,status,var,tmp
 call h5dclose_f(varid,status)

end subroutine

subroutine get_parameter_data_1d_real(grpid,var,nx,tmp)

 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx
 real,dimension(:),pointer :: tmp
 integer :: itile,varid,status
 integer(hsize_t) :: dims(1)
 real*8 :: tmp2(nx)
 allocate(tmp(nx))
 
 dims(1) = nx
 call h5dopen_f(grpid,var,varid,status)
 !call h5dread_f(varid,H5T_IEEE_F64LE,tmp,dims,status)
 call h5dread_f(varid,H5T_IEEE_F64LE,tmp2,dims,status)
 !call h5dread_f(varid,H5T_IEEE_F64LE,tmp2,dims,status)
 tmp = real(tmp2)
 !print*,status,var,tmp
 call h5dclose_f(varid,status) 

end subroutine 

subroutine get_parameter_data_2d_integer(grpid,var,nx,ny,tmp)

 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx,ny
 integer,dimension(:,:),pointer :: tmp
 integer :: itile,varid,status
 integer(hsize_t) :: dims(2)
 allocate(tmp(nx,ny))

 dims(1) = nx
 dims(2) = ny
 call h5dopen_f(grpid,var,varid,status)
 call h5dread_f(varid,H5T_STD_I32LE,tmp,dims,status)
 !print*,status,var,tmp
 call h5dclose_f(varid,status)

end subroutine

subroutine get_parameter_data_2d_real(grpid,var,nx,ny,tmp)

 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx,ny
 real,dimension(:,:),pointer :: tmp
 integer :: itile,varid,status
 integer(hsize_t) :: dims(2)
 real*8 :: tmp2(nx,ny)
 allocate(tmp(nx,ny))

 dims(1) = nx
 dims(2) = ny
 call h5dopen_f(grpid,var,varid,status)
 !call h5dread_f(varid,H5T_IEEE_F64LE,tmp,dims,status)
 call h5dread_f(varid,H5T_IEEE_F64LE,tmp2,dims,status)
 tmp = real(tmp2)
 !print*,status,var,tmp
 call h5dclose_f(varid,status)

end subroutine

end module predefined_tiles_mod
