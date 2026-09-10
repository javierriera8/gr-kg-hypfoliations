module numservicef90

  use prec

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First derivatives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdiff(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta2
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta2 = 1.0_wp/(2.0_wp*h)

if(.false.)then ! original in -r 452 - so far only implemented for second order of convergence
    if (bound.eq.'none') then
	first = 1
	last  = steps
     else if (bound.eq.'right') then
	first = 1
	last  = steps - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost
	last  = steps
    else if (bound.eq.'both') then
	first = 1 + nghost
	last  = steps - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
endif
if(.false..and.nghost.gt.1) then
    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - 1
    else if (bound.eq.'left') then
	first = 1 + nghost + 1
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + 1
	last  = steps - nghost - 1
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
    outvect(1+nghost-1)     =   (-3.0_wp*invect(1+nghost-1) + 4.0_wp*invect(2+nghost-1) - invect(3+nghost-1))*delta2
    outvect(steps-nghost+1) = -(-3.0_wp*invect(steps-nghost+1) + 4.0_wp*invect(steps-1-nghost+1) - invect(steps-2-nghost+1))*delta2
endif

    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost + nghost
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + nghost
	last  = steps - nghost - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if

!    outvect(1)     =   (-3.0_wp*invect(1) + 4.0_wp*invect(2) - invect(3))*delta2 ! is this correct?
!    outvect(steps) = -(-3.0_wp*invect(steps) + 4.0_wp*invect(steps-1) - invect(steps-2))*delta2

    outvect(1 + nghost)     =   (-3.0_wp*invect(1 + nghost) + 4.0_wp*invect(2 + nghost) - invect(3 + nghost))*delta2  ! is this correct? - let's see
    outvect(steps - nghost) = -(-3.0_wp*invect(steps - nghost) + 4.0_wp*invect(steps-1 - nghost) - invect(steps-2 - nghost))*delta2

    do i = first, last
!    do i = 2, steps - 1
       outvect(i) = ( invect(i+1) - invect(i-1) ) * delta2
    end do
  end subroutine vdiff


  subroutine mdiff(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call vdiff(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4th order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdiff_c4(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta12, delta2
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta12 = 1.0_wp/(12.0_wp*h)

if(.false.)then ! original in -r 510 - so far only implemented for second and fourth order of convergence
    if (bound.eq.'none') then
	first = 1
	last  = steps
     else if (bound.eq.'right') then
	first = 1
	last  = steps - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost
	last  = steps
    else if (bound.eq.'both') then
	first = 1 + nghost
	last  = steps - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
endif
    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost + nghost
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + nghost
	last  = steps - nghost - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if

!    outvect(1)     = (-25.0_wp*invect(1) + 48.0_wp*invect(2) - 36.0_wp*invect(3) + 16.0_wp*invect(4) - 3.0_wp*invect(5))*delta12
!    outvect(2)     = (- 3.0_wp*invect(1) - 10.0_wp*invect(2) + 18.0_wp*invect(3) -  6.0_wp*invect(4) +       invect(5))*delta12

!    outvect(steps-1)     = -(- 3.0_wp*invect(steps) - 10.0_wp*invect(steps-1) + 18.0_wp*invect(steps-2) -  6.0_wp*invect(steps-3) +       invect(steps-4))*delta12
!if(.true.) then
!    outvect(steps)       = -(-25.0_wp*invect(steps) + 48.0_wp*invect(steps-1) - 36.0_wp*invect(steps-2) + 16.0_wp*invect(steps-3) - 3.0_wp*invect(steps-4))*delta12
!else
!Decreasing order at outer boundary
!delta2 = 1.0_wp/(2.0_wp*h)
!outvect(steps) = -(-3.0_wp*invect(steps) + 4.0_wp*invect(steps-1)  - invect(steps-2))*delta2
!endif

    outvect(1 + nghost)     = (-25.0_wp*invect(1 + nghost) + 48.0_wp*invect(2 + nghost) - 36.0_wp*invect(3 + nghost) + 16.0_wp*invect(4 + nghost) - 3.0_wp*invect(5 + nghost))*delta12 ! is this correct? - let's see
    outvect(2 + nghost)     = (- 3.0_wp*invect(1 + nghost) - 10.0_wp*invect(2 + nghost) + 18.0_wp*invect(3 + nghost) -  6.0_wp*invect(4 + nghost) +       invect(5 + nghost))*delta12

    outvect(steps-1 - nghost)     = -(- 3.0_wp*invect(steps - nghost) - 10.0_wp*invect(steps-1 - nghost) + 18.0_wp*invect(steps-2 - nghost) -  6.0_wp*invect(steps-3 - nghost) +       invect(steps-4 - nghost))*delta12
if(.true.) then
    outvect(steps - nghost)       = -(-25.0_wp*invect(steps - nghost) + 48.0_wp*invect(steps-1 - nghost) - 36.0_wp*invect(steps-2 - nghost) + 16.0_wp*invect(steps-3 - nghost) - 3.0_wp*invect(steps-4 - nghost))*delta12
else
!Decreasing order at outer boundary
delta2 = 1.0_wp/(2.0_wp*h)
outvect(steps - nghost) = -(-3.0_wp*invect(steps - nghost) + 4.0_wp*invect(steps-1 - nghost)  - invect(steps-2 - nghost))*delta2
endif

!!!TEM
if(.true.)then
outvect(steps-1 - nghost) = ( -invect (-6 + steps - nghost) +7*invect (-5 + steps - nghost) -21*invect (-4 + steps - nghost) +36*invect (-3 + steps - nghost) - 43*invect (-2 + steps - nghost) +21*invect (-1 + steps - nghost) +1*invect (steps - nghost) )*delta12
outvect(steps - nghost) = ( invect (-6 + steps - nghost) -8*invect (-5 + steps - nghost) +28*invect (-4 + steps - nghost) -56*invect (-3 + steps - nghost) +71*invect (-2 + steps - nghost) -64*invect (-1 + steps - nghost) +28*invect (steps - nghost) )*delta12
endif

    do i = first, last
!    do i = 3, steps - 2
       outvect(i) = ( -invect(i+2) + 8.0_wp*invect(i+1) - 8.0_wp*invect(i-1) + invect(i-2) ) * delta12
    end do
  end subroutine vdiff_c4


  subroutine mdiff_c4(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call vdiff_c4(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiff_c4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 6th order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdiff_c6(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta60, delta2, delta12
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta60 = 1.0_wp/(60.0_wp*h)

if(.false.)then ! original in -r 511 - so far only implemented for second, fourth and sixth order of convergence
    if (bound.eq.'none') then
	first = 1
	last  = steps
     else if (bound.eq.'right') then
	first = 1
	last  = steps - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost
	last  = steps
    else if (bound.eq.'both') then
	first = 1 + nghost
	last  = steps - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
endif
    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost + nghost
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + nghost
	last  = steps - nghost - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if

!    outvect(1)     = (-147.0_wp*invect(1) + 360.0_wp*invect(2) - 450.0_wp*invect(3) + 400.0_wp*invect(4) - 225.0_wp*invect(5) + 72.0_wp*invect(6) - 10.0_wp*invect(7))*delta60
 !   outvect(2)     = (-10.0_wp*invect(1) - 77.0_wp*invect(2) + 150.0_wp*invect(3) - 100.0_wp*invect(4) + 50.0_wp*invect(5) - 15.0_wp*invect(6) + 2.0_wp*invect(7))*delta60
  !  outvect(3)     = (2.0_wp*invect(1) - 24.0_wp*invect(2) - 35.0_wp*invect(3) + 80.0_wp*invect(4) - 30.0_wp*invect(5) + 8.0_wp*invect(6) - invect(7))*delta60

   ! outvect(steps-2)=-(2.0_wp*invect(steps)-24.0_wp*invect(steps-1)-35.0_wp*invect(steps-2)+80.0_wp*invect(steps-3)-30.0_wp*invect(steps-4)+&
    !     & 8.0_wp*invect(steps-5)-invect(steps-6))*delta60
!if(.true.) then
 !   outvect(steps-1)=-(-10.0_wp*invect(steps)-77.0_wp*invect(steps-1)+150.0_wp*invect(steps-2)-100.0_wp*invect(steps-3)+50.0_wp*invect(steps-4)-&
  !       & 15.0_wp*invect(steps-5)+2.0_wp*invect(steps-6))*delta60
   ! outvect(steps)  = -(-147.0_wp*invect(steps) + 360.0_wp*invect(steps-1) - 450.0_wp*invect(steps-2) + 400.0_wp*invect(steps-3) -&
    !     & 225.0_wp*invect(steps-4) + 72.0_wp*invect(steps-5) - 10.0_wp*invect(steps-6))*delta60
!else
!Decreasing order at outer boundary
!delta12 = 1.0_wp/(12.0_wp*h)
!outvect(steps-1)     = -(- 3.0_wp*invect(steps) - 10.0_wp*invect(steps-1) + 18.0_wp*invect(steps-2) -  6.0_wp*invect(steps-3) +       invect(steps-4))*delta12
!delta2 = 1.0_wp/(2.0_wp*h)
!outvect(steps) = -(-3.0_wp*invect(steps) + 4.0_wp*invect(steps-1)  - invect(steps-2))*delta2
!endif

    outvect(1 + nghost)     = (-147.0_wp*invect(1 + nghost) + 360.0_wp*invect(2 + nghost) - 450.0_wp*invect(3 + nghost) + 400.0_wp*invect(4 + nghost) - 225.0_wp*invect(5 + nghost) + 72.0_wp*invect(6 + nghost) - 10.0_wp*invect(7 + nghost))*delta60 ! is this correct? - let's see
    outvect(2 + nghost)     = (-10.0_wp*invect(1 + nghost) - 77.0_wp*invect(2 + nghost) + 150.0_wp*invect(3 + nghost) - 100.0_wp*invect(4 + nghost) + 50.0_wp*invect(5 + nghost) - 15.0_wp*invect(6 + nghost) + 2.0_wp*invect(7 + nghost))*delta60
    outvect(3 + nghost)     = (2.0_wp*invect(1 + nghost) - 24.0_wp*invect(2 + nghost) - 35.0_wp*invect(3 + nghost) + 80.0_wp*invect(4 + nghost) - 30.0_wp*invect(5 + nghost) + 8.0_wp*invect(6 + nghost) - invect(7 + nghost))*delta60

    outvect(steps-2 - nghost)=-(2.0_wp*invect(steps - nghost)-24.0_wp*invect(steps-1 - nghost)-35.0_wp*invect(steps-2 - nghost)+80.0_wp*invect(steps-3 - nghost)-30.0_wp*invect(steps-4 - nghost) + 8.0_wp*invect(steps-5 - nghost)-invect(steps-6 - nghost))*delta60
    outvect(steps-1 - nghost)=-(-10.0_wp*invect(steps - nghost)-77.0_wp*invect(steps-1 - nghost)+150.0_wp*invect(steps-2 - nghost)-100.0_wp*invect(steps-3 - nghost)+50.0_wp*invect(steps-4 - nghost) - 15.0_wp*invect(steps-5 - nghost)+2.0_wp*invect(steps-6 - nghost))*delta60
    outvect(steps - nghost)  = -(-147.0_wp*invect(steps - nghost) + 360.0_wp*invect(steps-1 - nghost) - 450.0_wp*invect(steps-2 - nghost) + 400.0_wp*invect(steps-3 - nghost) - 225.0_wp*invect(steps-4 - nghost) + 72.0_wp*invect(steps-5 - nghost) - 10.0_wp*invect(steps-6 - nghost))*delta60

    do i = first, last
!    do i = 4, steps - 3
       outvect(i) = ( -invect(i-3) + invect(i+3) + 9.0_wp*(invect(i-2) - invect(i+2)) - 45.0_wp*(invect(i-1)-invect(i+1)) ) * delta60
    end do
  end subroutine vdiff_c6

  subroutine mdiff_c6(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call vdiff_c6(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiff_c6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 8th order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdiff_c8(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta840, delta2, delta12, delta60
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta840 = 1.0_wp/(840.0_wp*h)

if(.false.)then ! original in -r 519
    if (bound.eq.'none') then
	first = 1
	last  = steps
     else if (bound.eq.'right') then
	first = 1
	last  = steps - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost
	last  = steps
    else if (bound.eq.'both') then
	first = 1 + nghost
	last  = steps - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
endif
    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost + nghost
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + nghost
	last  = steps - nghost - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if

!    outvect(1)  =  ( -2283.0_wp*invect(1) + 6720.0_wp*invect(2) - 11760.0_wp*invect(3) + 15680.0_wp*invect(4) -  14700.0_wp*invect(5) +  9408.0_wp*invect(6) - 3920.0_wp*invect(7) + 960.0_wp*invect(8) - 105.0_wp*invect(9) )*delta840
 !   outvect(2)  = ( -105.0_wp*invect(1)-1338.0_wp*invect(2)+2940.0_wp*invect(3)-2940.0_wp*invect(4)+   2450.0_wp*invect(5)-1470.0_wp*invect(6)+   588.0_wp*invect(7) - 140.0_wp*invect(8) + 15.0_wp*invect(9) ) * delta840
  !  outvect(3)  = ( 15.0_wp*invect(1) - 240.0_wp*invect(2) - 798.0_wp*invect(3) + 1680.0_wp*invect(4) -  1050.0_wp*invect(5) + 560.0_wp*invect(6)-  210.0_wp*invect(7) + 48.0_wp*invect(8) - 5.0_wp*invect(9) ) * delta840
   ! outvect(4)  = (  -5.0_wp*invect(1) + 60.0_wp*invect(2) - 420.0_wp*invect(3) - 378.0_wp*invect(4) +   1050.0_wp*invect(5) - 420.0_wp*invect(6) +   140.0_wp*invect(7) - 30.0_wp*invect(8) + 3.0_wp*invect(9) ) * delta840

!    outvect(steps-3)  = - (  -5.0_wp*invect(steps) + 60.0_wp*invect(steps-1) - 420.0_wp*invect(steps-2) - 378.0_wp*invect(steps-3) + 1050.0_wp*invect(steps-4) - 420.0_wp*invect(steps-5) +  140.0_wp*invect(steps-6) - 30.0_wp*invect(steps-7) + 3.0_wp*invect(steps-8) ) * delta840
!if(.true.) then
 !   outvect(steps-2)  = - ( 15.0_wp*invect(steps) - 240.0_wp*invect(steps-1) - 798.0_wp*invect(steps-2) + 1680.0_wp*invect(steps-3) - 1050.0_wp*invect(steps-4) + 560.0_wp*invect(steps-5)-  210.0_wp*invect(steps-6) + 48.0_wp*invect(steps-7) - 5.0_wp*invect(steps-8) ) * delta840
  !  outvect(steps-1)  = - ( -105.0_wp*invect(steps)-1338.0_wp*invect(steps-1)+2940.0_wp*invect(steps-2)-2940.0_wp*invect(steps-3)+2450.0_wp*invect(steps-4)-1470.0_wp*invect(steps-5)+  588.0_wp*invect(steps-6) - 140.0_wp*invect(steps-7) + 15.0_wp*invect(steps-8) ) * delta840
   ! outvect(steps)  =  -( -2283.0_wp*invect(steps) + 6720.0_wp*invect(steps-1) - 11760.0_wp*invect(steps-2) + 15680.0_wp*invect(steps-3) - 14700.0_wp*invect(steps-4) +  9408.0_wp*invect(steps-5) - 3920.0_wp*invect(steps-6) + 960.0_wp*invect(steps-7) - 105.0_wp*invect(steps-8) )*delta840
!else
!Decreasing order at outer boundary
!delta60 = 1.0_wp/(60.0_wp*h)
!outvect(steps-2)=-(2.0_wp*invect(steps)-24.0_wp*invect(steps-1)-35.0_wp*invect(steps-2)+80.0_wp*invect(steps-3)-30.0_wp*invect(steps-4)+ 8.0_wp*invect(steps-5)-invect(steps-6))*delta60
!delta12 = 1.0_wp/(12.0_wp*h)
!outvect(steps-1)     = -(- 3.0_wp*invect(steps) - 10.0_wp*invect(steps-1) + 18.0_wp*invect(steps-2) -  6.0_wp*invect(steps-3) +       invect(steps-4))*delta12
!delta2 = 1.0_wp/(2.0_wp*h)
!outvect(steps) = -(-3.0_wp*invect(steps) + 4.0_wp*invect(steps-1)  - invect(steps-2))*delta2
!endif

    outvect(1 + nghost)  =  ( -2283.0_wp*invect(1 + nghost) + 6720.0_wp*invect(2 + nghost) - 11760.0_wp*invect(3 + nghost) + 15680.0_wp*invect(4 + nghost) -  14700.0_wp*invect(5 + nghost) +  9408.0_wp*invect(6 + nghost) - 3920.0_wp*invect(7 + nghost) + 960.0_wp*invect(8 + nghost) - 105.0_wp*invect(9 + nghost) )*delta840
    outvect(2 + nghost)  = ( -105.0_wp*invect(1 + nghost)-1338.0_wp*invect(2 + nghost)+2940.0_wp*invect(3 + nghost)-2940.0_wp*invect(4 + nghost)+   2450.0_wp*invect(5 + nghost)-1470.0_wp*invect(6 + nghost)+   588.0_wp*invect(7 + nghost) - 140.0_wp*invect(8 + nghost) + 15.0_wp*invect(9 + nghost) ) * delta840
    outvect(3 + nghost)  = ( 15.0_wp*invect(1 + nghost) - 240.0_wp*invect(2 + nghost) - 798.0_wp*invect(3 + nghost) + 1680.0_wp*invect(4 + nghost) -  1050.0_wp*invect(5 + nghost) + 560.0_wp*invect(6 + nghost)-  210.0_wp*invect(7 + nghost) + 48.0_wp*invect(8 + nghost) - 5.0_wp*invect(9 + nghost) ) * delta840
    outvect(4 + nghost)  = (  -5.0_wp*invect(1 + nghost) + 60.0_wp*invect(2 + nghost) - 420.0_wp*invect(3 + nghost) - 378.0_wp*invect(4 + nghost) +   1050.0_wp*invect(5 + nghost) - 420.0_wp*invect(6 + nghost) +   140.0_wp*invect(7 + nghost) - 30.0_wp*invect(8 + nghost) + 3.0_wp*invect(9 + nghost) ) * delta840

    outvect(steps-3 - nghost)  = - (  -5.0_wp*invect(steps - nghost) + 60.0_wp*invect(steps-1 - nghost) - 420.0_wp*invect(steps-2 - nghost) - 378.0_wp*invect(steps-3 - nghost) + 1050.0_wp*invect(steps-4 - nghost) - 420.0_wp*invect(steps-5 - nghost) +  140.0_wp*invect(steps-6 - nghost) - 30.0_wp*invect(steps-7 - nghost) + 3.0_wp*invect(steps-8 - nghost) ) * delta840
    outvect(steps-2 - nghost)  = - ( 15.0_wp*invect(steps - nghost) - 240.0_wp*invect(steps-1 - nghost) - 798.0_wp*invect(steps-2 - nghost) + 1680.0_wp*invect(steps-3 - nghost) - 1050.0_wp*invect(steps-4 - nghost) + 560.0_wp*invect(steps-5 - nghost)-  210.0_wp*invect(steps-6 - nghost) + 48.0_wp*invect(steps-7 - nghost) - 5.0_wp*invect(steps-8 - nghost) ) * delta840
    outvect(steps-1 - nghost)  = - ( -105.0_wp*invect(steps - nghost)-1338.0_wp*invect(steps-1 - nghost)+2940.0_wp*invect(steps-2 - nghost)-2940.0_wp*invect(steps-3 - nghost)+2450.0_wp*invect(steps-4 - nghost)-1470.0_wp*invect(steps-5 - nghost)+  588.0_wp*invect(steps-6 - nghost) - 140.0_wp*invect(steps-7 - nghost) + 15.0_wp*invect(steps-8 - nghost) ) * delta840
    outvect(steps - nghost)  =  -( -2283.0_wp*invect(steps - nghost) + 6720.0_wp*invect(steps-1 - nghost) - 11760.0_wp*invect(steps-2 - nghost) + 15680.0_wp*invect(steps-3 - nghost) - 14700.0_wp*invect(steps-4 - nghost) +  9408.0_wp*invect(steps-5 - nghost) - 3920.0_wp*invect(steps-6 - nghost) + 960.0_wp*invect(steps-7 - nghost) - 105.0_wp*invect(steps-8 - nghost) )*delta840

    do i = first, last
!    do i = 5, steps - 4
    outvect(i)=( 3.0_wp*(invect(i-4)-invect(4+i))+32.0_wp*(invect(3+i)-invect(i-3))+168.0_wp*(invect(i-2)-invect(2+i))+&
         & 672.0_wp*(invect(1 + i) - invect(i-1)) )*delta840
    end do
  end subroutine vdiff_c8

  subroutine mdiff_c8(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call vdiff_c8(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiff_c8


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Second derivatives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine v2diff(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta_2
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta_2 = 1.0_wp/(h*h)

    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost + nghost
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + nghost
	last  = steps - nghost - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if

if(.false..and.nghost.gt.1) then
    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - 1
    else if (bound.eq.'left') then
	first = 1 + nghost + 1
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + 1
	last  = steps - nghost - 1
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
    outvect(1+nghost-1)     =   (2*invect (1+nghost-1) - 5*invect (2+nghost-1) + 4*invect (3+nghost-1) - invect (4+nghost-1))*delta_2
    outvect(steps-nghost+1) = -((invect (-3 + steps-nghost+1) - 4*invect (-2 + steps-nghost+1) + 5*invect (-1 + steps-nghost+1) - 2*invect (steps-nghost+1))*delta_2)
endif

!    outvect(1)     = (2*invect (1) - 5*invect (2) + 4*invect (3) - invect (4))*delta_2
!    outvect(steps) = -((invect (-3 + steps) - 4*invect (-2 + steps) + 5*invect (-1 + steps) - 2*invect (steps))*delta_2)

    outvect(1 + nghost)     = (2*invect (1 + nghost) - 5*invect (2 + nghost) + 4*invect (3 + nghost) - invect (4 + nghost))*delta_2  ! is this correct? - let's see
    outvect(steps - nghost) = -((invect (-3 + steps - nghost) - 4*invect (-2 + steps - nghost) + 5*invect (-1 + steps - nghost) - 2*invect (steps - nghost))*delta_2)

    do i = first, last
!    do i = 2, steps - 1
       outvect(i) = (invect (-1 + i) - 2*invect (i) + invect (1 + i))*delta_2
    end do
  end subroutine v2diff


  subroutine m2diff(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call v2diff(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine m2diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4th order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine v2diff_c4(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta12_2, delta_2
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta12_2 = 1.0_wp/(12.0_wp*h*h)

if(.false.)then ! original in -r 510 - so far only implemented for second and fourth order of convergence
    if (bound.eq.'none') then
	first = 1
	last  = steps
     else if (bound.eq.'right') then
	first = 1
	last  = steps - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost
	last  = steps
    else if (bound.eq.'both') then
	first = 1 + nghost
	last  = steps - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
endif
    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost + nghost
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + nghost
	last  = steps - nghost - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if

!    outvect(1)     = (45*invect (1) - 154*invect (2) + 214*invect (3) - 156*invect (4) + 61*invect (5) - 10*invect (6))*delta12_2
!    outvect(2)     = (10*invect (1) - 15*invect (2) - 4*invect (3) + 14*invect (4) - 6*invect (5) + invect (6))*delta12_2

!    outvect(steps-1)     = (invect (-5 + steps) - 6*invect (-4 + steps) + 14*invect (-3 + steps) - 4*invect (-2 + steps) - 15*invect (-1 + steps) + 10*invect (steps))*delta12_2
!if(.true.) then
!    outvect(steps)       = (-10*invect (-5 + steps) + 61*invect (-4 + steps) - 156*invect (-3 + steps) + 214*invect (-2 + steps) - 154*invect (-1 + steps) + 45*invect (steps))*delta12_2
!else
!Decreasing order at outer boundary
!delta_2 = 1.0_wp/(h*h)
!outvect(steps) = -((invect (-3 + steps) - 4*invect (-2 + steps) + 5*invect (-1 + steps) - 2*invect (steps))*delta_2)
!endif

    outvect(1 + nghost)     = (45*invect (1 + nghost) - 154*invect (2 + nghost) + 214*invect (3 + nghost) - 156*invect (4 + nghost) + 61*invect (5 + nghost) - 10*invect (6 + nghost))*delta12_2 ! is this correct? - let's see
    outvect(2 + nghost)     = (10*invect (1 + nghost) - 15*invect (2 + nghost) - 4*invect (3 + nghost) + 14*invect (4 + nghost) - 6*invect (5 + nghost) + invect (6 + nghost))*delta12_2

    outvect(steps-1 - nghost)     = (invect (-5 + steps - nghost) - 6*invect (-4 + steps - nghost) + 14*invect (-3 + steps - nghost) - 4*invect (-2 + steps - nghost) - 15*invect (-1 + steps - nghost) + 10*invect (steps - nghost))*delta12_2
if(.true.) then
    outvect(steps - nghost)       = (-10*invect (-5 + steps - nghost) + 61*invect (-4 + steps - nghost) - 156*invect (-3 + steps - nghost) + 214*invect (-2 + steps - nghost) - 154*invect (-1 + steps - nghost) + 45*invect (steps - nghost))*delta12_2
else
!Decreasing order at outer boundary
delta_2 = 1.0_wp/(h*h)
outvect(steps - nghost) = -((invect (-3 + steps - nghost) - 4*invect (-2 + steps - nghost) + 5*invect (-1 + steps - nghost) - 2*invect (steps - nghost))*delta_2)
endif

!!!TEM
if(.true.)then
outvect(steps-1 - nghost) = ( -invect (-6 + steps - nghost) +7*invect (-5 + steps - nghost) -21*invect (-4 + steps - nghost) +34*invect (-3 + steps - nghost) - 19*invect (-2 + steps - nghost) -9*invect (-1 + steps - nghost) +9*invect (steps - nghost) )*delta12_2
outvect(steps - nghost) = ( 9*invect (-6 + steps - nghost) -64*invect (-5 + steps - nghost) +196*invect (-4 + steps - nghost) -336*invect (-3 + steps - nghost) +349*invect (-2 + steps - nghost) -208*invect (-1 + steps - nghost) +54*invect (steps - nghost) )*delta12_2
endif
!!!TEM


    do i = first, last
!    do i = 3, steps - 2
       outvect(i) = -(invect (-2 + i) - 16*invect (-1 + i) + 30*invect (i) - 16*invect (1 + i) + invect (2 + i))*delta12_2
    end do
  end subroutine v2diff_c4


  subroutine m2diff_c4(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call v2diff_c4(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine m2diff_c4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 6th order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine v2diff_c6(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta180_2, delta12_2, delta_2
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta180_2 = 1.0_wp/(180.0_wp*h*h)

if(.false.)then ! original in -r 511 - so far only implemented for second, fourth and sixth order of convergence
    if (bound.eq.'none') then
	first = 1
	last  = steps
     else if (bound.eq.'right') then
	first = 1
	last  = steps - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost
	last  = steps
    else if (bound.eq.'both') then
	first = 1 + nghost
	last  = steps - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
endif
    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost + nghost
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + nghost
	last  = steps - nghost - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if

!    outvect(1)     = (938*invect (1) - 4014*invect (2) + 7911*invect (3) - 9490*invect (4) + 7380*invect (5) - 3618*invect (6) + 1019*invect (7) - 126*invect (8))*delta180_2
 !   outvect(2)     = (126*invect (1) - 70*invect (2) - 486*invect (3) + 855*invect (4) - 670*invect (5) + 324*invect (6) - 90*invect (7) + 11*invect (8))*delta180_2
  !  outvect(3)     = (-11*invect (1) + 214*invect (2) - 378*invect (3) + 130*invect (4) + 85*invect (5) - 54*invect (6) + 16*invect (7) - 2*invect (8))*delta180_2

   ! outvect(steps-2)= (-2*invect (-7 + steps) + 16*invect (-6 + steps) - 54*invect (-5 + steps) + 85*invect (-4 + steps) + 130*invect (-3 + steps) - 378*invect (-2 + steps) + 214*invect (-1 + steps) - 11*invect (steps))*delta180_2
!if(.true.) then
 !   outvect(steps-1)= (11*invect (-7 + steps) - 90*invect (-6 + steps) + 324*invect (-5 + steps) - 670*invect (-4 + steps) + 855*invect (-3 + steps) - 486*invect (-2 + steps) - 70*invect (-1 + steps) + 126*invect (steps))*delta180_2
  !  outvect(steps)  = (-126*invect (-7 + steps) + 1019*invect (-6 + steps) - 3618*invect (-5 + steps) + 7380*invect (-4 + steps) - 9490*invect (-3 + steps) + 7911*invect (-2 + steps) - 4014*invect (-1 + steps) + 938*invect (steps))*delta180_2
!else
!Decreasing order at outer boundary
!delta12_2 = 1.0_wp/(12.0_wp*h*h)
!outvect(steps-1)     = (invect (-5 + steps) - 6*invect (-4 + steps) + 14*invect (-3 + steps) - 4*invect (-2 + steps) - 15*invect (-1 + steps) + 10*invect (steps))*delta12_2
!delta_2 = 1.0_wp/(h*h)
!outvect(steps) = -((invect (-3 + steps) - 4*invect (-2 + steps) + 5*invect (-1 + steps) - 2*invect (steps))*delta_2)
!endif

    outvect(1 + nghost)     = (938*invect (1 + nghost) - 4014*invect (2 + nghost) + 7911*invect (3 + nghost) - 9490*invect (4 + nghost) + 7380*invect (5 + nghost) - 3618*invect (6 + nghost) + 1019*invect (7 + nghost) - 126*invect (8 + nghost))*delta180_2
    outvect(2 + nghost)     = (126*invect (1 + nghost) - 70*invect (2 + nghost) - 486*invect (3 + nghost) + 855*invect (4 + nghost) - 670*invect (5 + nghost) + 324*invect (6 + nghost) - 90*invect (7 + nghost) + 11*invect (8 + nghost))*delta180_2
    outvect(3 + nghost)     = (-11*invect (1 + nghost) + 214*invect (2 + nghost) - 378*invect (3 + nghost) + 130*invect (4 + nghost) + 85*invect (5 + nghost) - 54*invect (6 + nghost) + 16*invect (7 + nghost) - 2*invect (8 + nghost))*delta180_2

    outvect(steps-2 - nghost)= (-2*invect (-7 + steps - nghost) + 16*invect (-6 + steps - nghost) - 54*invect (-5 + steps - nghost) + 85*invect (-4 + steps - nghost) + 130*invect (-3 + steps - nghost) - 378*invect (-2 + steps - nghost) + 214*invect (-1 + steps - nghost) - 11*invect (steps - nghost))*delta180_2
    outvect(steps-1 - nghost)= (11*invect (-7 + steps - nghost) - 90*invect (-6 + steps - nghost) + 324*invect (-5 + steps - nghost) - 670*invect (-4 + steps - nghost) + 855*invect (-3 + steps - nghost) - 486*invect (-2 + steps - nghost) - 70*invect (-1 + steps - nghost) + 126*invect (steps - nghost))*delta180_2
    outvect(steps - nghost)  = (-126*invect (-7 + steps - nghost) + 1019*invect (-6 + steps - nghost) - 3618*invect (-5 + steps - nghost) + 7380*invect (-4 + steps - nghost) - 9490*invect (-3 + steps - nghost) + 7911*invect (-2 + steps - nghost) - 4014*invect (-1 + steps - nghost) + 938*invect (steps - nghost))*delta180_2


    do i = first, last
!    do i = 4, steps - 3
       outvect(i) = (2*invect (-3 + i) - 27*invect (-2 + i) + 270*invect (-1 + i) - 490*invect (i) + 270*invect (1 + i) - 27*invect (2 + i) + 2*invect (3 + i))*delta180_2
    end do
  end subroutine v2diff_c6

  subroutine m2diff_c6(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call v2diff_c6(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine m2diff_c6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 8th order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine v2diff_c8(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta5040_2, delta_2, delta12_2, delta180_2
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta5040_2 = 1.0_wp/(5040.0_wp*h*h)

if(.false.)then ! original in -r 519
    if (bound.eq.'none') then
	first = 1
	last  = steps
     else if (bound.eq.'right') then
	first = 1
	last  = steps - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost
	last  = steps
    else if (bound.eq.'both') then
	first = 1 + nghost
	last  = steps - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
endif
    if (bound.eq.'none') then
	first = 1 + nghost
	last  = steps - nghost
     else if (bound.eq.'right') then
	first = 1 + nghost
	last  = steps - nghost - nghost
    else if (bound.eq.'left') then
	first = 1 + nghost + nghost
	last  = steps - nghost
    else if (bound.eq.'both') then
	first = 1 + nghost + nghost
	last  = steps - nghost - nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if

!    outvect(1)  =  (32575*invect (1) - 165924*invect (2) + 422568*invect (3) - 704368*invect (4) + 818874*invect (5) - 667800*invect (6) + 375704*invect (7) - 139248*invect (8) + 30663*invect (9) - 3044*invect (10))*delta5040_2
 !   outvect(2)  = (3044*invect (1) + 2135*invect (2) - 28944*invect (3) + 57288*invect (4) - 65128*invect (5) + 51786*invect (6) - 28560*invect (7) + 10424*invect (8) - 2268*invect (9) + 223*invect (10))*delta5040_2
  !  outvect(3)  = -(223*invect (1) - 5274*invect (2) + 7900*invect (3) + 2184*invect (4) - 10458*invect (5) + 8932*invect (6) - 4956*invect (7) + 1800*invect (8) - 389*invect (9) + 38*invect (10))*delta5040_2
   ! outvect(4)  = (38*invect (1) - 603*invect (2) + 6984*invect (3) - 12460*invect (4) + 5796*invect (5) + 882*invect (6) - 952*invect (7) + 396*invect (8) - 90*invect (9) + 9*invect (10))*delta5040_2

!    outvect(steps-3)  =  (9*invect (-9 + steps) - 90*invect (-8 + steps) + 396*invect (-7 + steps) - 952*invect (-6 + steps) + 882*invect (-5 + steps) + 5796*invect (-4 + steps) - 12460*invect (-3 + steps) + 6984*invect (-2 + steps) - 603*invect (-1 + steps) + 38*invect (steps))*delta5040_2
!if(.true.) then
 !   outvect(steps-2)  = -(38*invect (-9 + steps) - 389*invect (-8 + steps) + 1800*invect (-7 + steps) - 4956*invect (-6 + steps) + 8932*invect (-5 + steps) - 10458*invect (-4 + steps) + 2184*invect (-3 + steps) + 7900*invect (-2 + steps) - 5274*invect (-1 + steps) + 223*invect (steps))*delta5040_2
  !  outvect(steps-1)  = (223*invect (-9 + steps) - 2268*invect (-8 + steps) + 10424*invect (-7 + steps) - 28560*invect (-6 + steps) + 51786*invect (-5 + steps) - 65128*invect (-4 + steps) + 57288*invect (-3 + steps) - 28944*invect (-2 + steps) + 2135*invect (-1 + steps) + 3044*invect (steps))*delta5040_2
   ! outvect(steps)  = (-3044*invect (-9 + steps) + 30663*invect (-8 + steps) - 139248*invect (-7 + steps) + 375704*invect (-6 + steps) - 667800*invect (-5 + steps) + 818874*invect (-4 + steps) - 704368*invect (-3 + steps) + 422568*invect (-2 + steps) - 165924*invect (-1 + steps) + 32575*invect (steps))*delta5040_2
!else
!Decreasing order at outer boundary
!delta180_2 = 1.0_wp/(180.0_wp*h*h)
!outvect(steps-2)= (-2*invect (-7 + steps) + 16*invect (-6 + steps) - 54*invect (-5 + steps) + 85*invect (-4 + steps) + 130*invect (-3 + steps) - 378*invect (-2 + steps) + 214*invect (-1 + steps) - 11*invect (steps))*delta180_2
!delta12_2 = 1.0_wp/(12.0_wp*h*h)
!outvect(steps-1)     = (invect (-5 + steps) - 6*invect (-4 + steps) + 14*invect (-3 + steps) - 4*invect (-2 + steps) - 15*invect (-1 + steps) + 10*invect (steps))*delta12_2
!delta_2 = 1.0_wp/(h*h)
!outvect(steps) = -((invect (-3 + steps) - 4*invect (-2 + steps) + 5*invect (-1 + steps) - 2*invect (steps))*delta_2)
!endif

    outvect(1 + nghost)  =  (32575*invect (1 + nghost) - 165924*invect (2 + nghost) + 422568*invect (3 + nghost) - 704368*invect (4 + nghost) + 818874*invect (5 + nghost) - 667800*invect (6 + nghost) + 375704*invect (7 + nghost) - 139248*invect (8 + nghost) + 30663*invect (9 + nghost) - 3044*invect (10 + nghost))*delta5040_2
    outvect(2 + nghost)  = (3044*invect (1 + nghost) + 2135*invect (2 + nghost) - 28944*invect (3 + nghost) + 57288*invect (4 + nghost) - 65128*invect (5 + nghost) + 51786*invect (6 + nghost) - 28560*invect (7 + nghost) + 10424*invect (8 + nghost) - 2268*invect (9 + nghost) + 223*invect (10 + nghost))*delta5040_2
    outvect(3 + nghost)  = -(223*invect (1 + nghost) - 5274*invect (2 + nghost) + 7900*invect (3 + nghost) + 2184*invect (4 + nghost) - 10458*invect (5 + nghost) + 8932*invect (6 + nghost) - 4956*invect (7 + nghost) + 1800*invect (8 + nghost) - 389*invect (9 + nghost) + 38*invect (10 + nghost))*delta5040_2
    outvect(4 + nghost)  = (38*invect (1 + nghost) - 603*invect (2 + nghost) + 6984*invect (3 + nghost) - 12460*invect (4 + nghost) + 5796*invect (5 + nghost) + 882*invect (6 + nghost) - 952*invect (7 + nghost) + 396*invect (8 + nghost) - 90*invect (9 + nghost) + 9*invect (10 + nghost))*delta5040_2

    outvect(steps-3 - nghost)  =  (9*invect (-9 + steps - nghost) - 90*invect (-8 + steps - nghost) + 396*invect (-7 + steps - nghost) - 952*invect (-6 + steps - nghost) + 882*invect (-5 + steps - nghost) + 5796*invect (-4 + steps - nghost) - 12460*invect (-3 + steps - nghost) + 6984*invect (-2 + steps - nghost) - 603*invect (-1 + steps - nghost) + 38*invect (steps - nghost))*delta5040_2
    outvect(steps-2 - nghost)  = -(38*invect (-9 + steps - nghost) - 389*invect (-8 + steps - nghost) + 1800*invect (-7 + steps - nghost) - 4956*invect (-6 + steps - nghost) + 8932*invect (-5 + steps - nghost) - 10458*invect (-4 + steps - nghost) + 2184*invect (-3 + steps - nghost) + 7900*invect (-2 + steps - nghost) - 5274*invect (-1 + steps - nghost) + 223*invect (steps - nghost))*delta5040_2
    outvect(steps-1 - nghost)  = (223*invect (-9 + steps - nghost) - 2268*invect (-8 + steps - nghost) + 10424*invect (-7 + steps - nghost) - 28560*invect (-6 + steps - nghost) + 51786*invect (-5 + steps - nghost) - 65128*invect (-4 + steps - nghost) + 57288*invect (-3 + steps - nghost) - 28944*invect (-2 + steps - nghost) + 2135*invect (-1 + steps - nghost) + 3044*invect (steps - nghost))*delta5040_2
    outvect(steps - nghost)  = (-3044*invect (-9 + steps - nghost) + 30663*invect (-8 + steps - nghost) - 139248*invect (-7 + steps - nghost) + 375704*invect (-6 + steps - nghost) - 667800*invect (-5 + steps - nghost) + 818874*invect (-4 + steps - nghost) - 704368*invect (-3 + steps - nghost) + 422568*invect (-2 + steps - nghost) - 165924*invect (-1 + steps - nghost) + 32575*invect (steps - nghost))*delta5040_2

    do i = first, last
!    do i = 5, steps - 4
    outvect(i)= -(9*invect (-4 + i) - 128*invect (-3 + i) + 1008*invect (-2 + i) - 8064*invect (-1 + i) + 14350*invect (i) - 8064*invect (1 + i) + 1008*invect (2 + i) - 128*invect (3 + i) + 9*invect (4 + i))*delta5040_2
    end do
  end subroutine v2diff_c8

  subroutine m2diff_c8(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call v2diff_c8(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine m2diff_c8



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dissipation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdiff_diss4(invect, outvect, a, b, bound, nghost)

    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta16
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta16 = -1._wp/(16._wp*h) !changed sign


    ! Taken from vdiff2_c4, as nothing sensible needs to be done for the dissipation stuff
    outvect(1)       = 0
    outvect(2)       = ( invect(4) - 4d0*invect(3) + 6d0*invect(2) - 4d0*invect(1) + invect(3) ) * delta16 !0 ! take into account origin symmetry condition here even
    outvect(steps)   = 0
    outvect(steps-1) = ( 3*invect(steps)-3*invect(steps-1)+invect(steps-2) - 4d0*invect(steps) + 6d0*invect(steps-1) - 4d0*invect(steps-2) + invect(steps-3) ) * delta16 !0 ! seems to work

    do i = 3, steps - 2
       outvect(i) = ( invect(i+2) - 4d0*invect(i+1) + 6d0*invect(i) - 4d0*invect(i-1) + invect(i-2) ) * delta16
    end do
!    outvect(steps-2) = 0 !2*outvect(steps-3)-outvect(steps-4) ! special feature test
 !   outvect(steps-1) = 0 !2*outvect(steps-2)-outvect(steps-3) ! special feature test
!    outvect(steps-1) = (  5*invect(steps)-10*invect(steps-1)+10*invect(steps-2)-5*invect(steps-3)+invect(steps-4)      - 4d0*invect(steps) + 6d0*invect(steps-1) - 4d0*invect(steps-2) + invect(steps-3) ) * delta16 !test 2o
!    outvect(steps-4) = outvect(steps-4)*(3+0.5_wp)/(3+0.5_wp) ! special feature test
 !   outvect(steps-3) = outvect(steps-3)*(2+0.5_wp)/(3+0.5_wp) ! special feature test
  !  outvect(steps-2) = outvect(steps-2)*(1+0.5_wp)/(3+0.5_wp) ! special feature test
   ! outvect(steps-1) = outvect(steps-1)*(0+0.5_wp)/(3+0.5_wp) ! special feature test
  end subroutine vdiff_diss4


  subroutine vdiff_diss4odd(invect, outvect, a, b, bound, nghost) ! tried 04-06-2013
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    integer       :: steps, intervals
    real(kind=wp) :: h, delta16
    integer       :: i
    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta16 = -1._wp/(16._wp*h) !changed sign
    outvect(1)       = 0
    outvect(2)       = ( invect(4) - 4d0*invect(3) + 6d0*invect(2) - 4d0*invect(1) - invect(3) ) * delta16 !0 ! here odd
    outvect(steps)   = 0
    outvect(steps-1) = ( 3*invect(steps)-3*invect(steps-1)+invect(steps-2) - 4d0*invect(steps) + 6d0*invect(steps-1) - 4d0*invect(steps-2) + invect(steps-3) ) * delta16 !0 ! seems to work
    do i = 3, steps - 2
       outvect(i) = ( invect(i+2) - 4d0*invect(i+1) + 6d0*invect(i) - 4d0*invect(i-1) + invect(i-2) ) * delta16
    end do
!    outvect(steps-2) = 0 !2*outvect(steps-3)-outvect(steps-4) ! special feature test
 !   outvect(steps-1) = 0 !2*outvect(steps-2)-outvect(steps-3) ! special feature test
!    outvect(steps-1) = (  5*invect(steps)-10*invect(steps-1)+10*invect(steps-2)-5*invect(steps-3)+invect(steps-4)      - 4d0*invect(steps) + 6d0*invect(steps-1) - 4d0*invect(steps-2) + invect(steps-3) ) * delta16 !test 2o
!    outvect(steps-4) = outvect(steps-4)*(3+0.5_wp)/(3+0.5_wp) ! special feature test
 !   outvect(steps-3) = outvect(steps-3)*(2+0.5_wp)/(3+0.5_wp) ! special feature test
  !  outvect(steps-2) = outvect(steps-2)*(1+0.5_wp)/(3+0.5_wp) ! special feature test
   ! outvect(steps-1) = outvect(steps-1)*(0+0.5_wp)/(3+0.5_wp) ! special feature test
  end subroutine vdiff_diss4odd

  subroutine mdiff_diss4(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call vdiff_diss4(mat(j,:), tmp, a, b, bound, nghost)
       if(j.eq.6.or.j.eq.8.or.j.eq.9.or.j.eq.13) call vdiff_diss4odd(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiff_diss4


  subroutine vdiff_diss6(invect, outvect, a, b, bound, nghost)

    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta64
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta64 = 1._wp/(64._wp*h)

    outvect(1)       = 0
    outvect(2)       = 0
    outvect(3)       =  ( invect(6) - 6d0*invect(5) + 15d0*invect(4) - 20d0*invect(3) + 15d0*invect(2) -6d0*invect(1) + invect(5) ) * delta64 !0
    outvect(3) 		 = 0 ! for excision: also change next line
    outvect(3) 		 = ( invect(6) - 6d0*invect(5) + 15d0*invect(4) - 20d0*invect(3) + 15d0*invect(2) -6d0*invect(1) + (5*invect(1)-10*invect(2)+10*invect(3)-5*invect(4)+invect(5)) ) * delta64 ! test for extrapolation at the origin, basically following the expression for outvect(steps-2)
    outvect(steps)   = 0
    outvect(steps-1) = 0
!    outvect(steps-2) =  ( 3*invect(steps)-3*invect(steps-1)+invect(steps-2)    - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !0 ! first three terms copied form vdiff_diss4, seem to work
    outvect(steps-2) =  ( 5*invect(steps)-10*invect(steps-1)+10*invect(steps-2)-5*invect(steps-3)+invect(steps-4)      - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !0 ! better
!    outvect(steps-2) =  ( - invect(-5 + steps) + 6*invect(-4 + steps) - 15*invect(-3 + steps) + 20*invect(-2 + steps) - 15*invect(-1 + steps) + 6*invect(steps)  &
!&  - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !0
    do i = 4, steps - 3 !!-2
       outvect(i) = ( invect(i+3) - 6d0*invect(i+2) + 15d0*invect(i+1) - 20d0*invect(i) + 15d0*invect(i-1) -6d0*invect(i-2) + invect(i-3) ) * delta64
    end do
!outvect(steps-3)=0
!outvect(steps-4)=0
  !  outvect(steps-2) = ( 5*invect(steps)-10*invect(steps-1)+10*invect(steps-2)-5*invect(steps-3)+invect(steps-4)    - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !using the outflow condition for hD^5 to fill in the extra ghost point
!    outvect(1:10) = 0
!    outvect(steps-4) = 0 !2*outvect(steps-5)-outvect(steps-6) ! special feature test
 !   outvect(steps-3) = 0 !2*outvect(steps-4)-outvect(steps-5) ! special feature test
  !  outvect(steps-2) = 0 !2*outvect(steps-3)-outvect(steps-4) ! special feature test
!    outvect(steps-2) = ( invect(steps-6) - 7*invect(steps-5) + 21*invect(steps-4) - 35*invect(steps-3) + 35*invect(steps-2) - 21*invect(steps-1) + 7*invect(steps)   - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !test 2o ! 7th
!outvect(steps-2) = ( - invect(steps-5) + 6*invect(steps-4) - 15*invect(steps-3) + 20*invect(steps-2) - 15*invect(steps-1) + 6*invect(steps)    - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !test 2o ! 6th
!outvect(steps-2) = ( -invect(steps-7)+8*invect(steps-6) - 28*invect(steps-5)+56*invect(steps-4) - 70*invect(steps-3)+56*invect(steps-2) - 28*invect(steps-1)+8*invect(steps)    - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !test 2o ! 8th
!    outvect(steps-6) = outvect(steps-6)*(4+0.5_wp)/(4+0.5_wp) ! special feature test
 !   outvect(steps-5) = outvect(steps-5)*(3+0.5_wp)/(4+0.5_wp) ! special feature test
  !  outvect(steps-4) = outvect(steps-4)*(2+0.5_wp)/(4+0.5_wp) ! special feature test
   ! outvect(steps-3) = outvect(steps-3)*(1+0.5_wp)/(4+0.5_wp) ! special feature test
    !outvect(steps-2) = outvect(steps-2)*(0+0.5_wp)/(4+0.5_wp) ! special feature test
  end subroutine vdiff_diss6


  subroutine vdiff_diss6odd(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    integer       :: steps, intervals
    real(kind=wp) :: h, delta64
    integer       :: i
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta64 = 1._wp/(64._wp*h)
    outvect(1)       = 0
    outvect(2)       = 0
    outvect(3)       = ( invect(6) - 6d0*invect(5) + 15d0*invect(4) - 20d0*invect(3) + 15d0*invect(2) -6d0*invect(1) - invect(5) ) * delta64 !0
    outvect(3) 		 = 0 ! for excision: also change next line
    outvect(3) 		 = ( invect(6) - 6d0*invect(5) + 15d0*invect(4) - 20d0*invect(3) + 15d0*invect(2) -6d0*invect(1) + (5*invect(1)-10*invect(2)+10*invect(3)-5*invect(4)+invect(5)) ) * delta64 ! test for extrapolation at the origin, basically following the expression for outvect(steps-2)
    outvect(steps)   = 0
    outvect(steps-1) = 0
!    outvect(steps-2) =  ( 3*invect(steps)-3*invect(steps-1)+invect(steps-2)    - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !0 ! first three terms copied form vdiff_diss4, seem to work
    outvect(steps-2) =  ( 5*invect(steps)-10*invect(steps-1)+10*invect(steps-2)-5*invect(steps-3)+invect(steps-4)      - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !0 !
!    outvect(steps-2) =  ( - invect(-5 + steps) + 6*invect(-4 + steps) - 15*invect(-3 + steps) + 20*invect(-2 + steps) - 15*invect(-1 + steps) + 6*invect(steps)  &
!&  - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !0
    do i = 4, steps - 3 !!-2
       outvect(i) = ( invect(i+3) - 6d0*invect(i+2) + 15d0*invect(i+1) - 20d0*invect(i) + 15d0*invect(i-1) -6d0*invect(i-2) + invect(i-3) ) * delta64
    end do
!outvect(steps-3)=0
!outvect(steps-4)=0
  !  outvect(steps-2) = ( 5*invect(steps)-10*invect(steps-1)+10*invect(steps-2)-5*invect(steps-3)+invect(steps-4)    - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !using the outflow condition for hD^5 to fill in the extra ghost point
!    outvect(1:10) = 0
!    outvect(steps-4) = 0 !2*outvect(steps-5)-outvect(steps-6) ! special feature test
 !   outvect(steps-3) = 0 !2*outvect(steps-4)-outvect(steps-5) ! special feature test
  !  outvect(steps-2) = 0 !2*outvect(steps-3)-outvect(steps-4) ! special feature test
!    outvect(steps-2) = ( invect(steps-6) - 7*invect(steps-5) + 21*invect(steps-4) - 35*invect(steps-3) + 35*invect(steps-2) - 21*invect(steps-1) + 7*invect(steps)   - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !test 2o ! 7th
!outvect(steps-2) = ( - invect(steps-5) + 6*invect(steps-4) - 15*invect(steps-3) + 20*invect(steps-2) - 15*invect(steps-1) + 6*invect(steps)    - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !test 2o ! 6th
!outvect(steps-2) = ( -invect(steps-7)+8*invect(steps-6) - 28*invect(steps-5)+56*invect(steps-4) - 70*invect(steps-3)+56*invect(steps-2) - 28*invect(steps-1)+8*invect(steps)    - 6d0*invect(steps) + 15d0*invect(steps-1) - 20d0*invect(steps-2) + 15d0*invect(steps-3) -6d0*invect(steps-4) + invect(steps-5) ) * delta64 !test 2o ! 8th
!    outvect(steps-6) = outvect(steps-6)*(4+0.5_wp)/(4+0.5_wp) ! special feature test
 !   outvect(steps-5) = outvect(steps-5)*(3+0.5_wp)/(4+0.5_wp) ! special feature test
  !  outvect(steps-4) = outvect(steps-4)*(2+0.5_wp)/(4+0.5_wp) ! special feature test
   ! outvect(steps-3) = outvect(steps-3)*(1+0.5_wp)/(4+0.5_wp) ! special feature test
    !outvect(steps-2) = outvect(steps-2)*(0+0.5_wp)/(4+0.5_wp) ! special feature test
  end subroutine vdiff_diss6odd


  subroutine mdiff_diss6(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call vdiff_diss6(mat(j,:), tmp, a, b, bound, nghost)
       if(j.eq.6.or.j.eq.8.or.j.eq.9.or.j.eq.13) call vdiff_diss6odd(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiff_diss6


  subroutine vdiff_diss8(invect, outvect, a, b, bound, nghost)

    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta256
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta256 = -1.0_wp/(256._wp*h)

    outvect(1)       = 0
    outvect(2)       = 0
    outvect(3)       = 0
    outvect(4)       = ( invect(8) - 8._wp*invect(7) + 28._wp*invect(6) - 56._wp*invect(5) + 70._wp*invect(4) -56._wp*invect(3) + 28._wp*invect(2) - 8._wp*invect(1) + invect(7) )*delta256 !0
    outvect(steps)   = 0
    outvect(steps-1) = 0
    outvect(steps-2) = 0
    outvect(steps-3) =  0!( 3*invect(steps)-3*invect(steps-1)+invect(steps-2)      - 8._wp*invect(steps) + 28._wp*invect(steps-1) - 56._wp*invect(steps-2) + 70._wp*invect(steps-3) -56._wp*invect(steps-4) + 28._wp*invect(steps-5) - 8._wp*invect(steps-6) + invect(steps-7) )*delta256 !0 ! first three terms copied form vdiff_diss4

    do i = 5, steps - 4
       outvect(i) = ( invect(i+4) - 8._wp*invect(i+3) + 28._wp*invect(i+2) - 56._wp*invect(i+1) + 70._wp*invect(i) -56._wp*invect(i-1) + 28._wp*invect(i-2) - 8._wp*invect(i-3) + invect(i-4) )*delta256
    end do
  end subroutine vdiff_diss8

  subroutine vdiff_diss8odd(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    integer       :: steps, intervals
    real(kind=wp) :: h, delta256
    integer       :: i
    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta256 = -1.0_wp/(256._wp*h)
    outvect(1)       = 0
    outvect(2)       = 0
    outvect(3)       = 0
    outvect(4)       = ( invect(8) - 8._wp*invect(7) + 28._wp*invect(6) - 56._wp*invect(5) + 70._wp*invect(4) -56._wp*invect(3) + 28._wp*invect(2) - 8._wp*invect(1) - invect(7) )*delta256 !0
    outvect(steps)   = 0
    outvect(steps-1) = 0
    outvect(steps-2) = 0
    outvect(steps-3) =  0!( 3*invect(steps)-3*invect(steps-1)+invect(steps-2)      - 8._wp*invect(steps) + 28._wp*invect(steps-1) - 56._wp*invect(steps-2) + 70._wp*invect(steps-3) -56._wp*invect(steps-4) + 28._wp*invect(steps-5) - 8._wp*invect(steps-6) + invect(steps-7) )*delta256 !0 ! first three terms copied form vdiff_diss4
    do i = 5, steps - 4
       outvect(i) = ( invect(i+4) - 8._wp*invect(i+3) + 28._wp*invect(i+2) - 56._wp*invect(i+1) + 70._wp*invect(i) -56._wp*invect(i-1) + 28._wp*invect(i-2) - 8._wp*invect(i-3) + invect(i-4) )*delta256
    end do
  end subroutine vdiff_diss8odd

  subroutine mdiff_diss8(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call vdiff_diss8(mat(j,:), tmp, a, b, bound, nghost)
       if(j.eq.6.or.j.eq.8.or.j.eq.9.or.j.eq.13) call vdiff_diss8odd(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiff_diss8


  subroutine vdiff_diss10(invect, outvect, a, b, bound, nghost)

    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    integer       :: steps, intervals
    real(kind=wp) :: h, delta1024
    integer       :: i

    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta1024 = 1.0_wp/(1024._wp*h)

    outvect(1)       = 0
    outvect(2)       = 0
    outvect(3)       = 0
    outvect(4)       = 0
    outvect(5)       = ( (invect(10)) - 10.0_wp*(invect(9)+invect(1)) + 45.0_wp*(invect(8)+invect(2)) &
            & - 120.0_wp*(invect(7)+invect(3)) + 210.0_wp*(invect(6)+invect(4)) - 252.0_wp*invect(5) + invect(9))*delta1024 !0
    outvect(steps)   = 0
    outvect(steps-1) = 0
    outvect(steps-2) = 0
    outvect(steps-3) = 0
    outvect(steps-4) = 0

    do i = 6, steps - 5
       outvect(i) = ( (invect(i+5)+invect(i-5)) - 10.0_wp*(invect(i+4)+invect(i-4)) + 45.0_wp*(invect(i+3)+invect(i-3)) &
            & - 120.0_wp*(invect(i+2)+invect(i-2)) + 210.0_wp*(invect(i+1)+invect(i-1)) - 252.0_wp*invect(i) )*delta1024 !!?????????????????????????????????
    end do
  end subroutine vdiff_diss10

  subroutine vdiff_diss10odd(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    integer       :: steps, intervals
    real(kind=wp) :: h, delta1024
    integer       :: i
    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta1024 = 1.0_wp/(1024._wp*h)
    outvect(1)       = 0
    outvect(2)       = 0
    outvect(3)       = 0
    outvect(4)       = 0
    outvect(5)       = ( (invect(10)) - 10.0_wp*(invect(9)+invect(1)) + 45.0_wp*(invect(8)+invect(2)) &
            & - 120.0_wp*(invect(7)+invect(3)) + 210.0_wp*(invect(6)+invect(4)) - 252.0_wp*invect(5) - invect(9))*delta1024 !0
    outvect(steps)   = 0
    outvect(steps-1) = 0
    outvect(steps-2) = 0
    outvect(steps-3) = 0
    outvect(steps-4) = 0
    do i = 6, steps - 5
       outvect(i) = ( (invect(i+5)+invect(i-5)) - 10.0_wp*(invect(i+4)+invect(i-4)) + 45.0_wp*(invect(i+3)+invect(i-3)) &
            & - 120.0_wp*(invect(i+2)+invect(i-2)) + 210.0_wp*(invect(i+1)+invect(i-1)) - 252.0_wp*invect(i) )*delta1024 !!?????????????????????????????????
    end do
  end subroutine vdiff_diss10odd

  subroutine mdiff_diss10(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost

    integer	:: first, last

    real(kind=wp),                 dimension(size(mat, 2)) :: tmp

    integer       :: j

    do j = 1, size(mat, 1)
       call vdiff_diss10(mat(j,:), tmp, a, b, bound, nghost)
       if(j.eq.6.or.j.eq.8.or.j.eq.9.or.j.eq.13) call vdiff_diss10odd(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiff_diss10



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First one-sided derivatives (for advection terms)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdiffone(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    integer       :: steps, intervals
    real(kind=wp) :: h, delta2
    integer       :: i
    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta2 = 1.0_wp/(2.0_wp*h)

if(.false.)then
    if (bound.eq.'none'.or.bound.eq.'right') then
	first = 1 + nghost
    else if (bound.eq.'left'.or.bound.eq.'both') then
	first = 1 + nghost + nghost
    else
	stop 'The treatment of the stencil at the boundaries is not specified.'
    end if
	last  = steps - nghost
endif
	first = 1 + nghost + nghost
	last  = steps - nghost
!    outvect(1 + nghost)     =   (-3.0_wp*invect(1 + nghost) + 4.0_wp*invect(2 + nghost) - invect(3 + nghost))*delta2
        outvect(1 + nghost) = ( invect(1 + nghost+1) - invect(1 + nghost-1) ) * delta2
    do i = first, last
       outvect(i) = -(-3.0_wp*invect(i) + 4.0_wp*invect(i-1) - invect(i-2))*delta2
    end do
  end subroutine vdiffone


  subroutine mdiffone(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    real(kind=wp),                 dimension(size(mat, 2)) :: tmp
    integer       :: j
    do j = 1, size(mat, 1)
       call vdiffone(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiffone

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4th order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdiffone_c4(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    integer       :: steps, intervals
    real(kind=wp) :: h, delta12, delta2
    integer       :: i
    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta12 = 1.0_wp/(12.0_wp*h)

	first = 1 + nghost + nghost
	last  = steps - nghost
!    outvect(1 + nghost)     = (-25.0_wp*invect(1 + nghost) + 48.0_wp*invect(2 + nghost) - 36.0_wp*invect(3 + nghost) + 16.0_wp*invect(4 + nghost) - 3.0_wp*invect(5 + nghost))*delta12 ! original
!    outvect(2 + nghost)     = (- 3.0_wp*invect(1 + nghost) - 10.0_wp*invect(2 + nghost) + 18.0_wp*invect(3 + nghost) -  6.0_wp*invect(4 + nghost) +       invect(5 + nghost))*delta12 ! original
       outvect(1 + nghost) = ( -invect(1 + nghost+2) + 8.0_wp*invect(1 + nghost+1) - 8.0_wp*invect(1 + nghost-1) + invect(1 + nghost-2) ) * delta12 ! closest point to inner boundary centered
!       outvect(2 + nghost) = ( -invect(2 + nghost+2) + 8.0_wp*invect(2 + nghost+1) - 8.0_wp*invect(2 + nghost-1) + invect(2 + nghost-2) ) * delta12 !nope, this is centered and has to be offcentered
	outvect(steps - nghost)       = -(-25.0_wp*invect(steps - nghost) + 48.0_wp*invect(steps-1 - nghost) - 36.0_wp*invect(steps-2 - nghost) + 16.0_wp*invect(steps-3 - nghost) - 3.0_wp*invect(steps-4 - nghost))*delta12 !completely one-sided at the boundary
!    do i = first, last !+1 !completely one-sided at the boundary
    do i = first, last +1 !same one-point-on-the-right-offcentered like the rest at the boundary, use with centered stencil at the outer boundary
!       outvect(i)       = -(-25.0_wp*invect(i) + 48.0_wp*invect(i-1) - 36.0_wp*invect(i-2) + 16.0_wp*invect(i-3) - 3.0_wp*invect(i-4))*delta12 ! completely one-sided - have to keep one point on the right
        outvect(i-1)     = -(- 3.0_wp*invect(i) - 10.0_wp*invect(i-1) + 18.0_wp*invect(i-2) -  6.0_wp*invect(i-3) +  invect(i-4))*delta12 ! ok, keeping one point on the right
    end do
  end subroutine vdiffone_c4
  
    subroutine vdiffone_c4in(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    integer       :: steps, intervals
    real(kind=wp) :: h, delta12, delta2
    integer       :: i
    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta12 = 1.0_wp/(12.0_wp*h)

	first = 1 + nghost 
	last  = steps - nghost - nghost
	outvect(steps - nghost) = -( -invect(steps - nghost-2) + 8.0_wp*invect(steps - nghost-1) - 8.0_wp*invect(steps - nghost+1) + invect(steps - nghost+2) ) * delta12 ! closest point to outer boundary centered
	outvect(1 + nghost) = (-25.0_wp*invect(1 + nghost) + 48.0_wp*invect(1 + 1 + nghost) - 36.0_wp*invect(1 + 2 + nghost) + 16.0_wp*invect(1 + 3 + nghost) - 3.0_wp*invect(1 + 4 + nghost))*delta12 !completely one-sided at the inner boundary
    do i = first-1, last  !same one-point-on-the-right-offcentered like the rest at the boundary, use with centered stencil at the outer boundary
        outvect(i+1)     = (- 3.0_wp*invect(i) - 10.0_wp*invect(i+1) + 18.0_wp*invect(i+2) -  6.0_wp*invect(i+3) +  invect(i+4))*delta12 ! ok, keeping one point on the right
    end do
  end subroutine vdiffone_c4in

  subroutine mdiffone_c4(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    real(kind=wp),                 dimension(size(mat, 2)) :: tmp
    integer       :: j
    do j = 1, size(mat, 1)
       call vdiffone_c4(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiffone_c4
  
  subroutine mdiffone_c4in(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    real(kind=wp),                 dimension(size(mat, 2)) :: tmp
    integer       :: j
    do j = 1, size(mat, 1)
       call vdiffone_c4in(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiffone_c4in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 6th order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdiffone_c6(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    integer       :: steps, intervals
    real(kind=wp) :: h, delta60, delta2, delta12
    integer       :: i
    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta60 = 1.0_wp/(60.0_wp*h)

	first = 1 + nghost + nghost
	last  = steps - nghost
!    outvect(1 + nghost)     = (-147.0_wp*invect(1 + nghost) + 360.0_wp*invect(2 + nghost) - 450.0_wp*invect(3 + nghost) + 400.0_wp*invect(4 + nghost) - 225.0_wp*invect(5 + nghost) + 72.0_wp*invect(6 + nghost) - 10.0_wp*invect(7 + nghost))*delta60 ! is this correct? - let's see
!    outvect(2 + nghost)     = (-10.0_wp*invect(1 + nghost) - 77.0_wp*invect(2 + nghost) + 150.0_wp*invect(3 + nghost) - 100.0_wp*invect(4 + nghost) + 50.0_wp*invect(5 + nghost) - 15.0_wp*invect(6 + nghost) + 2.0_wp*invect(7 + nghost))*delta60
!    outvect(3 + nghost)     = (2.0_wp*invect(1 + nghost) - 24.0_wp*invect(2 + nghost) - 35.0_wp*invect(3 + nghost) + 80.0_wp*invect(4 + nghost) - 30.0_wp*invect(5 + nghost) + 8.0_wp*invect(6 + nghost) - invect(7 + nghost))*delta60
       outvect(1 + nghost) = ( -invect(1 + nghost-3) + invect(1 + nghost+3) + 9.0_wp*(invect(1 + nghost-2) - invect(1 + nghost+2)) - 45.0_wp*(invect(1 + nghost-1)-invect(1 + nghost+1)) ) * delta60
       outvect(2 + nghost) = ( -invect(2 + nghost-3) + invect(2 + nghost+3) + 9.0_wp*(invect(2 + nghost-2) - invect(2 + nghost+2)) - 45.0_wp*(invect(2 + nghost-1)-invect(2 + nghost+1)) ) * delta60
!       outvect(3 + nghost) = ( -invect(3 + nghost-3) + invect(3 + nghost+3) + 9.0_wp*(invect(3 + nghost-2) - invect(3 + nghost+2)) - 45.0_wp*(invect(3 + nghost-1)-invect(3 + nghost+1)) ) * delta60
    do i = first, last+1
!       outvect(i)  = -(-147.0_wp*invect(i) + 360.0_wp*invect(i-1) - 450.0_wp*invect(i-2) + 400.0_wp*invect(i-3) - 225.0_wp*invect(i-4) + 72.0_wp*invect(i-5) - 10.0_wp*invect(i-6))*delta60
    outvect(i-1)=-(-10.0_wp*invect(i)-77.0_wp*invect(i-1)+150.0_wp*invect(i-2)-100.0_wp*invect(i-3)+50.0_wp*invect(i-4) - 15.0_wp*invect(i-5)+2.0_wp*invect(i-6))*delta60
    end do
  end subroutine vdiffone_c6

  subroutine mdiffone_c6(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    real(kind=wp),                 dimension(size(mat, 2)) :: tmp
    integer       :: j
    do j = 1, size(mat, 1)
       call vdiffone_c6(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiffone_c6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 8th order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdiffone_c8(invect, outvect, a, b, bound, nghost)
    real(kind=wp), intent(in),  dimension(:)            :: invect
    real(kind=wp), intent(out), dimension(size(invect)) :: outvect
    real(kind=wp), intent(in)                           :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    integer       :: steps, intervals
    real(kind=wp) :: h, delta840, delta2, delta12, delta60
    integer       :: i
    ! executable statements
    steps = size(invect)
    intervals = steps-1
    h = (b-a)/dble(intervals)
    delta840 = 1.0_wp/(840.0_wp*h)

	first = 1 + nghost + nghost
	last  = steps - nghost

!    outvect(1 + nghost)  =  ( -2283.0_wp*invect(1 + nghost) + 6720.0_wp*invect(2 + nghost) - 11760.0_wp*invect(3 + nghost) + 15680.0_wp*invect(4 + nghost) -  14700.0_wp*invect(5 + nghost) +  9408.0_wp*invect(6 + nghost) - 3920.0_wp*invect(7 + nghost) + 960.0_wp*invect(8 + nghost) - 105.0_wp*invect(9 + nghost) )*delta840
!    outvect(2 + nghost)  = ( -105.0_wp*invect(1 + nghost)-1338.0_wp*invect(2 + nghost)+2940.0_wp*invect(3 + nghost)-2940.0_wp*invect(4 + nghost)+   2450.0_wp*invect(5 + nghost)-1470.0_wp*invect(6 + nghost)+   588.0_wp*invect(7 + nghost) - 140.0_wp*invect(8 + nghost) + 15.0_wp*invect(9 + nghost) ) * delta840
!    outvect(3 + nghost)  = ( 15.0_wp*invect(1 + nghost) - 240.0_wp*invect(2 + nghost) - 798.0_wp*invect(3 + nghost) + 1680.0_wp*invect(4 + nghost) -  1050.0_wp*invect(5 + nghost) + 560.0_wp*invect(6 + nghost)-  210.0_wp*invect(7 + nghost) + 48.0_wp*invect(8 + nghost) - 5.0_wp*invect(9 + nghost) ) * delta840
!    outvect(4 + nghost)  = (  -5.0_wp*invect(1 + nghost) + 60.0_wp*invect(2 + nghost) - 420.0_wp*invect(3 + nghost) - 378.0_wp*invect(4 + nghost) +   1050.0_wp*invect(5 + nghost) - 420.0_wp*invect(6 + nghost) +   140.0_wp*invect(7 + nghost) - 30.0_wp*invect(8 + nghost) + 3.0_wp*invect(9 + nghost) ) * delta840
    outvect(1 + nghost)=( 3.0_wp*(invect(1 + nghost-4)-invect(4+1 + nghost))+32.0_wp*(invect(3+1 + nghost)-invect(1 + nghost-3))+168.0_wp*(invect(1 + nghost-2)-invect(2+1 + nghost))+672.0_wp*(invect(1 + 1 + nghost) - invect(1 + nghost-1)) )*delta840
    outvect(2 + nghost)=( 3.0_wp*(invect(2 + nghost-4)-invect(4+2 + nghost))+32.0_wp*(invect(3+2 + nghost)-invect(2 + nghost-3))+168.0_wp*(invect(2 + nghost-2)-invect(2+2 + nghost))+672.0_wp*(invect(1 + 2 + nghost) - invect(2 + nghost-1)) )*delta840
    outvect(3 + nghost)=( 3.0_wp*(invect(3 + nghost-4)-invect(4+3 + nghost))+32.0_wp*(invect(3+3 + nghost)-invect(3 + nghost-3))+168.0_wp*(invect(3 + nghost-2)-invect(2+3 + nghost))+672.0_wp*(invect(1 + 3 + nghost) - invect(3 + nghost-1)) )*delta840
!    outvect(4 + nghost)=( 3.0_wp*(invect(4 + nghost-4)-invect(4+4 + nghost))+32.0_wp*(invect(3+4 + nghost)-invect(4 + nghost-3))+168.0_wp*(invect(4 + nghost-2)-invect(2+4 + nghost))+672.0_wp*(invect(1 + 4 + nghost) - invect(4 + nghost-1)) )*delta840

    do i = first, last+1
!    outvect(i)  =  -( -2283.0_wp*invect(i) + 6720.0_wp*invect(i-1) - 11760.0_wp*invect(i-2) + 15680.0_wp*invect(i-3) - 14700.0_wp*invect(i-4) +  9408.0_wp*invect(i-5) - 3920.0_wp*invect(i-6) + 960.0_wp*invect(i-7) - 105.0_wp*invect(i-8) )*delta840
    outvect(i-1)  = - ( -105.0_wp*invect(i)-1338.0_wp*invect(i-1)+2940.0_wp*invect(i-2)-2940.0_wp*invect(i-3)+2450.0_wp*invect(i-4)-1470.0_wp*invect(i-5)+  588.0_wp*invect(i-6) - 140.0_wp*invect(i-7) + 15.0_wp*invect(i-8) ) * delta840
    end do
  end subroutine vdiffone_c8

  subroutine mdiffone_c8(mat, a, b, bound, nghost)
    real(kind=wp), intent(inout),  dimension(:,:)          :: mat
    real(kind=wp), intent(in)                              :: a, b
     character(len = 42):: bound
     integer 		:: nghost
    integer	:: first, last
    real(kind=wp),                 dimension(size(mat, 2)) :: tmp
    integer       :: j
    do j = 1, size(mat, 1)
       call vdiffone_c8(mat(j,:), tmp, a, b, bound, nghost)
       mat(j,:) = tmp
    end do
  end subroutine mdiffone_c8


end module numservicef90
