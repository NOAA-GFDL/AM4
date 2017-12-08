module tile_diag_buff_mod

implicit none
private

! ==== public interfaces =====================================================
public :: diag_buff_type
public :: new_diag_buff, delete_diag_buff, realloc_diag_buff
! ==== end of public interfaces ==============================================
interface new_diag_buff
   module procedure diag_buff_ctor
   module procedure diag_buff_copy_ctor
end interface

! storage for tile diagnostic data
type :: diag_buff_type
   real   , pointer :: data(:) => NULL()
   logical, pointer :: mask(:) => NULL()
end type diag_buff_type

! ==== module constants =====================================================
integer, parameter :: MIN_DIAG_BUFF_SIZE = 1

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
function diag_buff_ctor() result(buffer)
  type(diag_buff_type), pointer :: buffer

  integer :: m ! initial size of the buffer

  allocate(buffer)
  m = MIN_DIAG_BUFF_SIZE
  allocate(buffer%mask(m),buffer%data(m))
  ! initialize buffer content
  buffer%mask(:) = .FALSE.
  buffer%data(:) = 0.0
end function


! ============================================================================
function diag_buff_copy_ctor(buffer) result(ptr)
  type(diag_buff_type), pointer :: ptr ! return value
  type(diag_buff_type), intent(in) :: buffer ! buffer to copy

  allocate(ptr)
  allocate(ptr%mask(size(buffer%mask)),ptr%data(size(buffer%data)))
  ! initialize buffer content
  ptr%mask(:) = buffer%mask(:)
  ptr%data(:) = buffer%data(:)
end function


! ============================================================================
subroutine delete_diag_buff(buffer)
  type(diag_buff_type), pointer :: buffer

  if(.not.associated(buffer)) return
  deallocate(buffer%mask,buffer%data)
  deallocate(buffer)

end subroutine


! ============================================================================
! reallocates buffer to have at least m elements
subroutine realloc_diag_buff(buffer, m)
  type(diag_buff_type), intent(inout) :: buffer
  integer             , intent(in)    :: m

  real    , pointer :: new_data(:)
  logical , pointer :: new_mask(:)
  integer           :: n

  ! n is size of the original buffer; m is the current size of the buffer
  ! for all diagnostic fields
  n = size(buffer%data(:))
  ! do nothing if buffer is big enough
  if(n >= m) return

  allocate(new_data(m), new_mask(m))
  new_data(1:n) = buffer%data(1:n) ; new_data(n+1:m) = 0.0
  new_mask(1:n) = buffer%mask(1:n) ; new_mask(n+1:m) = .FALSE.
  deallocate(buffer%data, buffer%mask)

  buffer%data=>new_data
  buffer%mask=>new_mask

end subroutine

! =============================================================================
subroutine merge_diag_buffs(t1,w1,t2,w2)
  type(diag_buff_type), intent(in)    :: t1
  type(diag_buff_type), intent(inout) :: t2
  real, intent(in) :: w1, w2 ! relative weights

  ! ---- local vars
  real :: x1, x2 ! normalized relative weights
  real :: y1, y2
  integer :: i

  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  do i = 1,size(t2%data)
     y1 = 0.0 ; if (t1%mask(i)) y1 = t1%data(i)
     y2 = 0.0 ; if (t2%mask(i)) y2 = t2%data(i)
     t2%data(i) = y1*x1 + y2*x2
     t2%mask(i) = t1%mask(i).or.t1%mask(i)
  enddo
end subroutine

end module tile_diag_buff_mod
