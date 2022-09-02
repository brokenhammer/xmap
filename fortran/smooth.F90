!###### MODULE SMOOTH_FUNCTION ###########################
! Contains generic functions for smoothing

module smooth_function
implicit none

public
contains

! Generic 1D smooth function
! Proxy to different types of 1D smoothing
! Arguments:
!     array: the data needs smoothing, will be smoothen IN PLACE
!     n: length of array
!     boundary_condition: optional, default to 0
!       chosen boundary condition:
!              0 for normal fix boundary, 
!              1 for periodic boundary.
!     smooth_type: optional default to 0
!       smoothing algorithm for choose:
!              0 for "121" smoothing,
!              1 for "14641" smoothing.

subroutine smooth_1d(array, n, boundary_condition, smooth_type)

  real(8), intent(inout):: array(:)
  integer, intent(in):: n
  integer, intent(in), optional:: boundary_condition
  integer, intent(in), optional:: smooth_type
  integer :: bc, st ! local copies of the boundary_condition and 
                    ! smooth_type arguments 
 
  if (present(smooth_type)) then
    st = smooth_type
  else
    st = 0 ! default to use 1-2-1 smoothing function
  endif

  if (present(boundary_condition)) then
    bc = boundary_condition
  else
    bc = 0 ! default to use normal fix boundary
  endif

! Now switch through the smoothing types and boundary conditions

  select case (st)
    case (0)
      if (bc==0) then
        call smooth_fix_121(array, n)
      else if (bc==1) then
        call smooth_periodic_121(array, n)
      else
        print *, "Invalid BC option in smooth_1d: ", bc
        stop
      endif
    case (1)
      if (bc==0) then
        call smooth_fix_14641(array, n)
      else if (bc==1) then
        call smooth_periodic_14641(array, n)
      else
        print *, "Invalid BC option in smooth_1d: ", bc
        stop
      endif
    case default
      print *,"Invalid smoothing type in smooth_1d: ", st
      stop
  end select
 
end subroutine smooth_1d

! subroutine to smooth the fix boundary data using 121 formula
subroutine smooth_fix_121(array, n)
  real(8), intent(inout) :: array(:)
  integer, intent (in) :: n
  integer :: i
  real(8) :: temp(n)
  temp(:) = array(:)
  do i=2, n-1
     array(i) = 0.25_8 *(temp(i-1) + 2.0_8*temp(i) + temp(i+1))
  enddo
end subroutine smooth_fix_121

! subroutine to smooth the periodic data using 121 formula
subroutine smooth_periodic_121(array, n)
  real(8), intent(inout) :: array(:)
  integer, intent (in) :: n
  integer :: i
  real(8) :: temp(n)
! if input array is not periodic, the array(n) will be assigned by array(1)
  if (ABS(array(1)-array(n))> SPACING(array(1))+SPACING(array(n)) ) then ! test if array(1) equals array(n)
    write(*,*) "SMOOTH_1D WARNING: periodic boundary condition &
	&applied to non-periodic data."//NEW_LINE('a')//"array(n)=array(1) is enforced."
    array(n)=array(1)
  endif
  temp(:) = array(:)
  do i=2, n-1
     array(i) = 0.25_8 *(temp(i-1) + 2.0_8*temp(i) + temp(i+1))
  enddo
! smooth end points using periodic boundary condition
  array(1) = 0.25_8 *(temp(n-1) + 2.0_8 *temp(1) + temp(2))
  array(n) = array(1)
end subroutine smooth_periodic_121

! subroutine to smooth the fix boundary data using 14641 formula
subroutine smooth_fix_14641(array, n)
  real(8), intent(inout) :: array(:)
  integer, intent (in) :: n
  integer :: i
  real(8) :: temp(n), coeff
  temp(:) = array(:)
  coeff = 1.0_8/12.0_8
  do i=3, n-2
     array(i) = coeff * (-(temp(i-2)+temp(i+2)) + 4.0_8*(temp(i-1)+temp(i+1)) + 6.0_8*temp(i))
  enddo
end subroutine smooth_fix_14641

! subroutine to smooth the periodic data using 14641 formula
subroutine smooth_periodic_14641(array, n)
  real(8), intent(inout) :: array(:)
  integer, intent (in) :: n
  integer :: i
  real(8) :: temp(n), coeff
! if input array is not periodic, the array(n) will be assigned by array(1)
  if (ABS(array(1)-array(n))> SPACING(array(1))+SPACING(array(n)) ) then ! test if array(1) equals array(n)
    write(*,*) "SMOOTH_1D WARNING: periodic boundary condition &
	&applied to non-periodic data."//NEW_LINE('a')//"array(n)=array(1) is enforced."
    array(n)=array(1)
  endif 
  temp(:) = array(:)
  coeff = 1.0_8/12.0_8
  do i=3, n-2
     array(i) = coeff * (-(temp(i-2)+temp(i+2)) + 4.0_8*(temp(i-1)+temp(i+1)) + 6.0_8*temp(i))
  enddo
! smooth end points using periodic boundary condition
  array(1) = coeff * (-(temp(n-2)+temp(3)) + 4.0_8*(temp(n-1)+temp(2)) + 6.0_8*temp(1))
  array(n) = array(1)
  array(2) = coeff * (-(temp(n-1)+temp(4)) + 4.0_8*(temp(1)+temp(3)) + 6.0_8*temp(2))
  array(n-1) = coeff * (-(temp(n-3)+temp(2)) + 4.0_8*(temp(n-2)+temp(1)) + 6.0_8*temp(n-1))
end subroutine smooth_periodic_14641

end module smooth_function

