! -*-f90-*-
! $Id$

! some sanity checks
#ifndef F90_TYPE
#error F90_TYPE is not defined: must be one of FORTRAN 90 types
#endif

#ifndef NF_TYPE
#error NF_TYPE is not defined: must be netcdf type name corresponding to F90_TYPE
#endif

#ifndef READ_0D_FPTR
#error name of subroutine READ_0D_FPTR is not defined
#endif

#ifndef WRITE_0D_FPTR
#error name of subroutine WRITE_0D_FPTR is not defined
#endif

! ============================================================================
subroutine READ_0D_FPTR(ncid,name,fptr,rec)
  integer           , intent(in) :: ncid ! netcdf id
  character(len=*)  , intent(in) :: name ! name of the variable to read
  integer, optional , intent(in) :: rec  ! record number (in case there are 
                                         ! several in the file) 
  ! subroutine returning the pointer to the data to be written
  interface ; subroutine fptr(cohort, ptr)
     use vegn_cohort_mod, only : vegn_cohort_type
     type(vegn_cohort_type), pointer :: cohort ! input
     F90_TYPE, pointer :: ptr ! returned pointer to the data
   end subroutine fptr
  end interface

  ! ---- local constants
  character(*), parameter :: module_name = 'read_cohort_data_r0d_fptr'

  ! ---- local vars
  integer :: i, n
  integer :: rec_     ! record number
  integer :: ntiles   ! size of the tile dimension in restart file
  integer :: ncohorts ! total number of cohorts in restart file
  integer :: bufsize  ! size of the input buffer
  integer :: idxid    ! id of the index dimension
  integer :: start(1),count(1) ! definition of slab for reading
  integer, allocatable :: idx(:) ! index dimension
  F90_TYPE, allocatable :: data(:) ! data to be read
  F90_TYPE, pointer :: ptr ! pointer to the individual cohort data
  type(vegn_cohort_type), pointer :: cohort

  ! assign the internal record number
  if(present(rec)) then
     rec_ = rec
  else
     rec_ = 1
  endif

  ! get the size of the tile dimension
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! get the length of cohort compressed index
  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))
  __NF_ASRT__(nfu_inq_var(ncid,cohort_index_name,id=idxid))

  ! allocate data
  bufsize=min(input_buf_size,ncohorts)
  allocate(data(bufsize),idx(bufsize))

  do n = 1, ncohorts, bufsize
     ! read the cohort index
     __NF_ASRT__(nf_get_vara_int(ncid,idxid,n,min(bufsize,ncohorts-n+1),idx))
     ! read the data
     start(1) = n; count(1) = min(bufsize,ncohorts-n+1)
     __NF_ASRT__(nfu_get_rec(ncid,name,rec_,data,start,count))
     
     ! distribute data over cohorts
     do i = 1, size(idx)
        call get_cohort_by_idx ( idx(i), lnd%nlon, lnd%nlat, ntiles,&
             land_tile_map, lnd%ls, lnd%gs, lnd%ge, cohort)
        if (associated(cohort)) then
           call fptr(cohort, ptr)
           if(associated(ptr)) ptr = data(i)
        endif
     enddo
  enddo
  
  ! free allocated memory
  deallocate(data,idx)
  
end subroutine READ_0D_FPTR

! ============================================================================
subroutine WRITE_0D(ncid,name,data,long_name,units,record)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  F90_TYPE        , intent(in) :: data(:) ! data to be written
  character(len=*), intent(in), optional :: units, long_name
  integer         , intent(in), optional :: record

  ! ---- local vars
  integer :: iret, varid, record_, p, k
  F90_TYPE, allocatable :: buffer(:) ! input buffer for data from other PEs
  integer,  allocatable :: nc(:) ! number of cohorts per PE in IO domain
  integer :: dimids(2), ndims

  ! if this processor isn't the root IO processor, simply send data to the root 
  ! IO processor and return from the subroutine
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(size(data), plen=1,          to_pe=lnd%io_pelist(1), tag=COMM_TAG_1)
     call mpp_send(data(1),    plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_2)
  else
     allocate(nc(size(lnd%io_pelist)))
     nc(1) = size(data)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(nc(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_1)
     enddo
     ! gather data from the processors in io_domain
     allocate(buffer(sum(nc(:))))
     buffer(1:nc(1)) = data(:)
     k=nc(1)+1
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buffer(k), glen=nc(p), from_pe=lnd%io_pelist(p), tag=COMM_TAG_2)
        k = k+nc(p)
     enddo

     ! create variable, if it does not exist
     if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
        ! get the ID of cohort dimension
        __NF_ASRT__(nf_inq_dimid(ncid,cohort_index_name,dimids(1)))
        
        ndims = 1
        if(present(record)) then
           if(nf_inq_unlimdim(ncid,dimids(2))==NF_NOERR) then
              ndims = 2
           endif
        endif
        __NF_ASRT__(nfu_def_var(ncid,name,NF_TYPE,dimids(1:ndims),long_name,units))
     endif
     ! write data
     iret = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
     record_ = 1
     if(present(record)) record_ = record
     __NF_ASRT__(nfu_put_rec(ncid,name,record_,buffer))
     deallocate(buffer)
  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to 
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
  ! free allocated memory
end subroutine WRITE_0D