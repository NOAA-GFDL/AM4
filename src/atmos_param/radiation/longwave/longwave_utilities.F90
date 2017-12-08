 
module longwave_utilities_mod


use fms_mod,      only: write_version_number, &
                        error_mesg, FATAL
!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!-------  version number --------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------

public :: longwave_utilities_init, &
          longwave_utilities_end, &
          locate_in_table,  &
          looktab, table_alloc

interface looktab
    module procedure  looktab_type1, looktab_type2, looktab_type3
end interface

interface table_alloc
   module procedure    table1_alloc, table2_alloc, table3_alloc
end interface

!---------------------------------------------------------------------
!------- public derived-types ------

public :: gas_tf_type

type gas_tf_type
     real, dimension(:,:,:),   pointer :: tdav=>NULL(),    &
                                          tlsqu=>NULL(),   &
                                          tmpdiff=>NULL(), &
                                          tstdav=>NULL(),  &
                                          n2o9c=>NULL(),   &
                                          tn2o17=>NULL()
     real, dimension(:,:,:),   pointer :: co2nbl=>NULL()
     real, dimension(:,:,:),   pointer :: co2990nbl=>NULL(), &
                                          co2900nbl=>NULL(), &
                                          co21070nbl=>NULL()
     real, dimension(:,:,:,:), pointer :: co2spnb=>NULL()
     real, dimension(:,:,:),   pointer :: co2990spnb=>NULL()
     real, dimension(:,:,:),   pointer :: co2900spnb=>NULL()
     real, dimension(:,:,:),   pointer :: co21070spnb=>NULL()
     real, dimension(:,:),     pointer :: a1=>NULL(),    &
                                          a2=>NULL()
end type gas_tf_type

!------------------------------------------------------------------

public longwave_tables1_type

type longwave_tables1_type
    real, dimension(:,:), pointer  ::  vae=>NULL(),   &
                                       td=>NULL(), &
                                       md=>NULL(), &
                                       cd=>NULL()
end type longwave_tables1_type

!--------------------------------------------------------------------

public longwave_tables2_type

type longwave_tables2_type
    real, dimension(:,:,:), pointer  ::  vae=>NULL(),  &
                                         td=>NULL(),  &
                                         md=>NULL(),   &
                                         cd=>NULL()
end type longwave_tables2_type

!---------------------------------------------------------------------

public longwave_tables3_type

type longwave_tables3_type
     real,  dimension(:,:), pointer    ::  vae=>NULL(),   &
                                           td=>NULL()
end type longwave_tables3_type

!---------------------------------------------------------------------

public lw_clouds_type

type lw_clouds_type
     real, dimension(:,:,:,:),   pointer :: taucld_rndlw=>NULL(), &
                                            taucld_mxolw=>NULL(), &
                                            taunbl_mxolw=>NULL()
end type lw_clouds_type

!------------------------------------------------------------------

public lw_table_type

type lw_table_type
     real, dimension(:),    pointer :: bdlocm=>NULL(),   &
                                       bdhicm=>NULL(),  &
                                       bandlo=>NULL(),  &
                                       bandhi=>NULL()
     integer, dimension(:), pointer :: iband=>NULL()
end type lw_table_type

!------------------------------------------------------------------

public optical_path_type

type optical_path_type
     real, dimension (:,:,:,:), pointer :: empl1f=>NULL(),  &
                                           empl2f=>NULL(),  &
                                           vrpfh2o=>NULL(), &
                                           xch2obd=>NULL(),  &
                                           tphfh2o=>NULL(), &
                                           avephif=>NULL(), &
                                           totaerooptdep=>NULL()
     real, dimension (:,:,:),   pointer :: empl1=>NULL(), &
                                           empl2=>NULL(),  &
                                           var1=>NULL(), &
                                           var2=>NULL(), &
                                           emx1f=>NULL(),   &
                                           emx2f=>NULL(),   &
                                           totvo2=>NULL(),  &
                                           avephi=>NULL(),&
                                           totch2obdwd=>NULL(), &
                                           xch2obdwd=>NULL(), &
                                           totphi=>NULL(),   &
                                           cntval=>NULL(), &
                                           toto3=>NULL(),   &
                                           tphio3=>NULL(),  &
                                           var3=>NULL(),  &
                                           var4=>NULL(),        &
                                           wk=>NULL(),         &
                                           rh2os=>NULL(),  &
                                           rfrgn=>NULL(),  &
                                           tfac=>NULL(), &
                                           totaerooptdep_15=>NULL(), &
                                           totf11=>NULL(),   &
                                           totf12=>NULL(),  &
                                           totf113=>NULL(),   &
                                           totf22=>NULL()
      real, dimension (:,:), pointer    :: emx1=>NULL(),  &
                                           emx2=>NULL(),  &
                                           csfah2o=>NULL(), &
                                           aerooptdep_KE_15=>NULL()
end type optical_path_type

!------------------------------------------------------------------

public sealw99_control_type

type sealw99_control_type
    character(len=16) :: continuum_form
    character(len=16) :: linecatalog_form
    logical           :: do_ch4lbltmpint
    logical           :: do_n2olbltmpint
    logical           :: do_lwcldemiss
    logical           :: do_h2o
    logical           :: do_o3
    logical           :: do_ch4
    logical           :: do_n2o
    logical           :: do_co2
    logical           :: do_co2_10um
    logical           :: do_cfc
end type sealw99_control_type

!------------------------------------------------------------------

public table_axis_type

type table_axis_type
  integer :: first_col
  real    :: min_val
  real    :: max_val
  real    :: tab_inc
end type table_axis_type

!---------------------------------------------------------------------
!------- public data -------

type (table_axis_type),        public   ::    &
               temp_1 = table_axis_type( 1, 100.0, 370.0, 10.0), &
               mass_1 = table_axis_type( 1, -16.0,   1.9,  0.1)

type (sealw99_control_type),  public   :: &
      Sealw99_control = sealw99_control_type(  '    ', '    ', &
                                              .false., .false., .false., &
                                              .false., .false., .false., &
                                              .false., .false., .false., &
                                              .false. )

!---------------------------------------------------------------------
!------- private data -------

logical :: module_is_initialized=.false.   ! module is initialized ?

!---------------------------------------------------------------------

CONTAINS

!#####################################################################
!
!                     PUBLIC SUBROUTINES
!
!#####################################################################

subroutine longwave_utilities_init

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    write version number to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)

!-------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------

end subroutine longwave_utilities_init

!#####################################################################

subroutine longwave_utilities_end

!--------------------------------------------------------------------
!    this is the destructor for longwave_utilities_mod
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilites_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!-------------------------------------------------------------------

end subroutine longwave_utilities_end

!#####################################################################
! <SUBROUTINE NAME="locate_in_table">
!  <OVERVIEW>
!   Subroutine to locate index and residual value from an array provided 
!   with array and axis information
!  </OVERVIEW>
!  <DESCRIPTION>
!     given array x and an arithmetic sequence of table column headings
!     tabxmin, tabxmin+tabdeltax, ..., corresponding to column ixlow, 
!     ixlow+1, ..., ixupp, Locate returns the array ix is column 
!     indices and the array dx of residuals.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call locate_in_table(table_axis, x, dx, ix, k_min, k_max)
!  </TEMPLATE>
!  <IN NAME="table_axis" TYPE="table_axis_type">
!   table_axis contains the axis information such as, min, increment,
!   and first column values.
!  </IN>
!  <IN NAME="x" TYPE="real">
!   array from which data is to be searched
!  </IN>
!  <OUT NAME="dx" TYPE="real">
!   residual between x and x(ix+first_column)
!  </OUT>
!  <OUT NAME="ix" TYPE="integer">
!   index values of the searched domain in the array
!  </OUT>
!  <IN NAME="k_min" TYPE="integer">
!   minimum k value of the search domain 
!  </IN>
!  <IN NAME="k_max" TYPE="integer">
!   maximum k value of the search domain
!  </IN>
! </SUBROUTINE>

subroutine locate_in_table (table_axis, x, dx, ix, k_min, k_max)

!---------------------------------------------------------------------
!    given array x and an arithmetic sequence of table column headings
!    tabxmin, tabxmin+tabdeltax, ..., corresponding to column ixlow, 
!    ixlow+1, ..., ixupp, locate_in_table returns the array ix of
!    column indices and the array dx of residuals.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!----------------------------------------------------------------------

type(table_axis_type),     intent(in)  :: table_axis
real,    dimension(:,:,:), intent(in)  :: x
integer,                   intent(in)  :: k_min, k_max
real,    dimension(:,:,:), intent(out) :: dx
integer, dimension(:,:,:), intent(out) :: ix

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    table_axis
!    x
!    k_min
!    k_max
!
!  intent(out) variables:
!
!    dx
!    ix
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(x,1), size(x,2), size(x,3))  ::  fx
      integer     ::  k

!---------------------------------------------------------------------
!  local variables:
!
!     fx
!     table_min
!     table_inc
!     k
!     table_col
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do k=k_min,k_max
        fx (:,:,k) = AINT((x(:,:,k) - table_axis%min_val )/  &
                     table_axis%tab_inc)
        ix (:,:,k) = INT(fx(:,:,k)) + table_axis%first_col
        dx (:,:,k) = x(:,:,k) - fx(:,:,k)*table_axis%tab_inc - &
                     table_axis%min_val
      end do

!---------------------------------------------------------------------

end subroutine locate_in_table

!####################################################################
! <SUBROUTINE NAME="looktab_type1">
!  <OVERVIEW>
!   Subroutine to calculate answer from input differentials.
!  </OVERVIEW>
!  <DESCRIPTION>
!   given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!     arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!     y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!     from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call looktab_type1 (tab, ix, iy, dx, dy, answer, k_min, k_max)
!  </TEMPLATE>
!  <IN NAME="tab" TYPE="longwave_tables1_type">
!   The data array that contains function values and differentials
!  </IN>
!  <IN NAME="ix" TYPE="integer">
!   x subscript of input data array
!  </IN>
!  <IN NAME="iy" TYPE="integer">
!   y subscript of input data array
!  </IN>
!  <IN NAME="dx" TYPE="real">
!   x step in the x subscript space
!  </IN>
!  <IN NAME="dy" TYPE="real">
!   y step in the y subscript space
!  </IN>
!  <OUT NAME="answer" TYPE="real">
!   the answer to be calculated
!  </OUT>
!  <IN NAME="k_min" TYPE="integer">
!   the minimum k value of the domain
!  </IN>
!  <IN NAME="k_max" TYPE="integer">
!   the maximum k value of the domain
!  </IN>
! </SUBROUTINE>

subroutine looktab_type1 (tab, ix, iy, dx, dy, answer, k_min, k_max)

!----------------------------------------------------------------------
!    given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!    arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!    y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!    from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!--------------------------------------------------------------------

type(longwave_tables1_type), intent(in)  :: tab
integer,dimension(:,:,:),    intent(in)  :: ix, iy
real,   dimension(:,:,:),    intent(in)  :: dx, dy
real,   dimension(:,:,:),    intent(out) :: answer
integer,                     intent(in)  :: k_min, k_max

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    tab
!    ix
!    iy
!    dx
!    dy
!    k_min
!    k_max
!
!  intent(out) variables:
!
!    answer
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer    ::  i_min, i_max, j_min, j_max, i, j, k

!--------------------------------------------------------------------
!  local variables:
!
!    i_min
!    i_max
!    j_min
!    j_max
!    i,j,k
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
        do j=j_min, j_max
          do i=i_min, i_max
            answer(i,j,k) =                                         &
                                      tab%vae (ix(i,j,k), iy(i,j,k)) + &
                            dx(i,j,k)*tab%td  (ix(i,j,k), iy(i,j,k)) + &
                            dy(i,j,k)*tab%md  (ix(i,j,k), iy(i,j,k)) + &
                  dx(i,j,k)*dy(i,j,k)*tab%cd(ix(i,j,k), iy(i,j,k))
          end do
        end do
      end do

!---------------------------------------------------------------------

end subroutine looktab_type1

!#####################################################################
! <SUBROUTINE NAME="looktab_type2">
!  <OVERVIEW>
!   Subroutine to calculate answer from input differentials.
!  </OVERVIEW>
!  <DESCRIPTION>
!   given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!     arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!     y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!     from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!     The difference between this version about the version above is
!     that the differential arrays are 3 dimensional.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call looktab_type2 (tab, ix, iy, dx, dy, answer, k_min, k_max, m)
!  </TEMPLATE>
!  <IN NAME="tab" TYPE="longwave_tables2_type">
!   The data array that contains function values and differentials
!  </IN>
!  <IN NAME="ix" TYPE="integer">
!   x subscript of input data array
!  </IN>
!  <IN NAME="iy" TYPE="integer">
!   y subscript of input data array
!  </IN>
!  <IN NAME="dx" TYPE="real">
!   x step in the x subscript space
!  </IN>
!  <IN NAME="dy" TYPE="real">
!   y step in the y subscript space
!  </IN>
!  <OUT NAME="answer" TYPE="real">
!   the answer to be calculated
!  </OUT>
!  <IN NAME="k_min" TYPE="integer">
!   the minimum k value of the domain
!  </IN>
!  <IN NAME="k_max" TYPE="integer">
!   the maximum k value of the domain
!  </IN>
!  <IN NAME="m" TYPE="integer">
!   the z indice of the differential arrays
!  </IN>
! </SUBROUTINE>

subroutine looktab_type2 (tab, ix, iy, dx, dy, answer, k_min, k_max, m)

!-------------------------------------------------------------------
!    given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!    arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!    y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!    from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!--------------------------------------------------------------------

type(longwave_tables2_type), intent(in)   :: tab
integer, dimension (:,:,:),  intent(in)   :: ix, iy
integer,                     intent(in)   :: m
real, dimension (:,:,:),     intent(in)   :: dx, dy
real, dimension (:,:,:),     intent(out)  :: answer
integer,                     intent(in)   :: k_min, k_max

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    tab
!    ix
!    iy
!    m
!    dx
!    dy
!    k_min
!    k_max
!
!  intent(out) variables:
!
!    answer
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

       integer    ::    i_min, i_max, j_min, j_max
       integer    ::    i, j, k

!--------------------------------------------------------------------
!  local variables:
!
!    i_min
!    i_max
!    j_min
!    j_max
!    i,j,k
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
        do j=j_min, j_max
          do i=i_min, i_max
            answer(i,j,k) =                                           &
                                   tab%vae (ix(i,j,k), iy(i,j,k),m) + &
                         dx(i,j,k)*tab%td (ix(i,j,k), iy(i,j,k),m) + &
                         dy(i,j,k)*tab%md (ix(i,j,k), iy(i,j,k),m) + &
               dx(i,j,k)*dy(i,j,k)*tab%cd   (ix(i,j,k), iy(i,j,k),m)
           end do
        end do
      end do

!--------------------------------------------------------------------

end subroutine looktab_type2

!###################################################################
! <SUBROUTINE NAME="looktab_type3">
!  <OVERVIEW>
!   Subroutine to calculate answer from input differentials.
!  </OVERVIEW>
!  <DESCRIPTION>
!   given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!     arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!     y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!     from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!   In this version, only f(x,y) and f(x,y)+dx*df/dx is used. Probably
!   the f(x,y) is homogeneous in y space.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call looktab_type3 (tab, ix, dx,  answer, k_min, k_max, n)
!  </TEMPLATE>
!  <IN NAME="tab" TYPE="longwave_tables3_type">
!   The data array that contains function values and differentials
!  </IN>
!  <IN NAME="ix" TYPE="integer">
!   x subscript of input data array
!  </IN>
!  <IN NAME="dx" TYPE="real">
!   x step in the x subscript space
!  </IN>
!  <OUT NAME="answer" TYPE="real">
!   the answer to be calculated
!  </OUT>
!  <IN NAME="k_min" TYPE="integer">
!   the minimum k value of the domain
!  </IN>
!  <IN NAME="k_max" TYPE="integer">
!   the maximum k value of the domain
!  </IN>
!  <IN NAME="n" TYPE="integer">
!   the z indice of the differential arrays
!  </IN>
! </SUBROUTINE>

subroutine looktab_type3 (tab, ix, dx,  answer, k_min, k_max, n)

!----------------------------------------------------------------------
!
!    given arrays ix(:,:,:) and dx(:,:,:) of integer subscripts and!
!    differences from x(:,:,:) and constant column subscript iyconst, 
!    calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:)) from four tables
!    of values f, df/dx, df/dy, and d2f/dxdy.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!-----------------------------------------------------------------------

type(longwave_tables3_type), intent(in)  :: tab
integer, dimension (:,:,:),  intent(in)  :: ix
integer,                     intent(in)  :: n
real,    dimension(:,:,:),   intent(in)  :: dx
real,    dimension(:,:,:),   intent(out) :: answer
integer,                     intent(in)  :: k_min, k_max

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    tab
!    ix
!    n
!    dx
!    k_min
!    k_max
!
!  intent(out) variables:
!
!    answer
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer    :: i_min, i_max, j_min, j_max
      integer    :: i, j, k

!--------------------------------------------------------------------
!  local variables:
!
!    i_min
!    i_max
!    j_min
!    j_max
!    i,j,k
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!-----------------------------------------------------------------
      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
        do j=j_min, j_max
          do i=i_min, i_max
                answer(i,j,k) =                                 &
                                      tab%vae (ix(i,j,k),n) +   &
                            dx(i,j,k)*tab%td(ix(i,j,k),n)
          end do
        end do
      end do

!------------------------------------------------------------------

end subroutine  looktab_type3

!#####################################################################
! <SUBROUTINE NAME="table1_alloc">
!  <OVERVIEW>
!   Allocate the longwave tables.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Allocate the longwave tables based on 2 dimension sizes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call table1_alloc(tab, dim1, dim2)
!  </TEMPLATE>
!  <INOUT NAME="tab" TYPE="longwave_tables1_type">
!   The longwave tables
!  </INOUT>
!  <IN NAME="dim1" TYPE="integer">
!   size of the x dimension
!  </IN>
!  <IN NAME="dim2" TYPE="integer">
!   size of the y dimension
!  </IN>
! </SUBROUTINE>

subroutine table1_alloc (tab, dim1, dim2)

!------------------------------------------------------------------
!    table1_alloc allocates the arrays contained in a 
!    longwave_tables1_type variable.
!------------------------------------------------------------------

type(longwave_tables1_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2

!-------------------------------------------------------------------
! intent(in) variables:
!
!     dim1
!     dim2
!
!  intent(inout) variables:
!
!     tab
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate the component arrays of a longwave_tables1_type variable.
!---------------------------------------------------------------------
      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))
      allocate (tab%md (dim1, dim2))
      allocate (tab%cd (dim1, dim2))

!---------------------------------------------------------------------

end subroutine table1_alloc

!####################################################################
! <SUBROUTINE NAME="table2_alloc">
!  <OVERVIEW>
!   Allocate the longwave tables.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Allocate the longwave tables based on 3 dimension sizes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call table2_alloc(tab, dim1, dim2)
!  </TEMPLATE>
!  <INOUT NAME="tab" TYPE="longwave_tables2_type">
!   The longwave tables
!  </INOUT>
!  <IN NAME="dim1" TYPE="integer">
!   size of the x dimension
!  </IN>
!  <IN NAME="dim2" TYPE="integer">
!   size of the y dimension
!  </IN>
!  <IN NAME="dim3" TYPE="integer">
!   size of the z dimension
!  </IN>
! </SUBROUTINE>

subroutine table2_alloc (tab, dim1, dim2, dim3)

!------------------------------------------------------------------
!    table2_alloc allocates the arrays contained in a 
!    longwave_tables2_type variable.
!------------------------------------------------------------------

type(longwave_tables2_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2, dim3

!-------------------------------------------------------------------
! intent(in) variables:
!
!     dim1
!     dim2
!     dim3
!
!  intent(inout) variables:
!
!     tab
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate the component arrays of a longwave_tables2_type variable.
!---------------------------------------------------------------------
      allocate (tab%vae(dim1, dim2, dim3))
      allocate (tab%td (dim1, dim2, dim3))
      allocate (tab%md (dim1, dim2, dim3))
      allocate (tab%cd (dim1, dim2, dim3))

!--------------------------------------------------------------------

end subroutine table2_alloc


!#####################################################################
! <SUBROUTINE NAME="table3_alloc">
!  <OVERVIEW>
!   Allocate the longwave tables.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Allocate the longwave tables based on 2 dimension sizes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call table3_alloc(tab, dim1, dim2)
!  </TEMPLATE>
!  <INOUT NAME="tab" TYPE="longwave_tables3_type">
!   The longwave tables
!  </INOUT>
!  <IN NAME="dim1" TYPE="integer">
!   size of the x dimension
!  </IN>
!  <IN NAME="dim2" TYPE="integer">
!   size of the y dimension
!  </IN>
! </SUBROUTINE>

subroutine table3_alloc (tab, dim1, dim2)

type(longwave_tables3_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2

!-------------------------------------------------------------------
! intent(in) variables:
!
!     dim1
!     dim2
!
!  intent(inout) variables:
!
!     tab
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate the component arrays of a longwave_tables3_type variable.
!---------------------------------------------------------------------
      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))

end subroutine table3_alloc

!##################################################################

end module longwave_utilities_mod

