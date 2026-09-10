module output

  use prec
  use parameters_mod
  use molparameters

  implicit none

contains

 subroutine outputparam_testing(parameters)
    implicit none

    type(parameters_type), intent(inout) :: parameters

    if(parameters%output%every_t_step <= 0) then
       parameters%output%every_t_step = 1
    endif

    if(parameters%output%every_x_step <= 0) then
       parameters%output%every_x_step = 1
    endif
  end subroutine outputparam_testing

  subroutine set_default_output_parameters
    implicit none

    default_parameters%output%every_t_step = 1
    default_parameters%output%every_x_step = 1
    default_parameters%output%suffix       = '.yg' !for pygraph
  end subroutine set_default_output_parameters


  subroutine write_out(time, index, xval, yval, filename, nghost, &
 &                     sample_factor)
    implicit none

    ! All arguments are declared with intent IN to avoid accidental overwrites.
    real(kind=wp),                     intent(in) :: time
    integer,                           intent(in) :: index, nghost
    real(kind=wp), dimension(:),       intent(in) :: xval, yval
    character(len=*),                  intent(in) :: filename
    integer,                           intent(in) :: sample_factor

    logical, save               :: first = .true.
    integer, parameter          :: unit = 3
    integer, save               :: dim, first_index

    integer ::  count, first_point, last_point

    first_point = lbound(xval, 1) + nghost
    last_point  = ubound(xval, 1) - nghost

    if (first) then
      first_index = index

      if (first_point == last_point) then
         dim = 0 ! do 0D-output only
         write(*,*) 'problem is 0-dimensional, switching to 0D-output'
      else
         dim = 1 ! do 1D-output
         write(*,*) 'problem is 1-dimensional, switching to 1D-output'
      endif
    endif

     if (index == first_index) then
         open(unit, file = trim(filename), status = 'replace', err = 120)
     else
         open(unit, file = trim(filename), status = 'unknown', err = 120, &
    &         position = 'append')
     end if

    if (dim == 0) then  ! 0D - output
       write(unit, *) time, yval
    else                ! 1D - output
       write(unit, *)
       write(unit, *)

       ! Write time information and then the actual (optionally downsampled)
       ! data. Add two empty lines to indicate the end of the time slice.
       write(unit,'(a8,f16.10)') '"Time = ', time

       do count = first_point, last_point, sample_factor
          write(unit,'(2(e37.28e3))') xval(count), yval(count)
       end do

    endif
    close (unit) ! Close output file.

    first = .false.
    return

120 call file_error(filename)
  end subroutine write_out


  subroutine write_out_par(k, index, xval, yval, filename, parameters)!!!!!!!!!!!!!!1111

    use parameters_mod

    implicit none

    ! All arguments are declared with intent IN to avoid accidental overwrites.

    integer,                           intent(in) :: index, k !!!!!!!!!!1
    real(kind=wp), dimension(:),       intent(in) :: xval, yval
    character(len=*),                  intent(in) :: filename
    character(len=15) :: tim !for pygraph

    type(parameters_type),             intent(in) :: parameters

    integer :: stagdiff(1:9)

    character(len=128)          :: fullname
    logical, save               :: first = .true.
    integer, parameter          :: unit = 3
    integer, save               :: dim, first_index

    integer ::  count, first_point, last_point, coeflp

    if (parameters%grid%origin.eq.'nostag') then
	stagdiff = (/0,0,0,0,0,0,0,0,0/)
	coeflp = 0
    elseif (parameters%grid%origin.eq.'stag') then
	stagdiff = (/0,1,4,13,40,121,364,1093,3280/) ! to select the corresponding points (staggered grid at x = 0)
	coeflp = 1
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!stagdiff=0 ! set to zero to output the closest point to the left boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    fullname = trim(file_prefix) // trim(filename) // trim(parameters%output%suffix)

! Change here to output ghost points
    first_point = lbound(xval, 1) + parameters%grid%nghost
    last_point  = ubound(xval, 1) - parameters%grid%nghost

    if (first) then
      first_index = index

      if (first_point == last_point) then
         dim = 0 ! do 0D-output only
         write(*,*) 'problem is 0-dimensional, switching to 0D-output'
      else
         dim = 1 ! do 1D-output
         write(*,*) 'problem is 1-dimensional, switching to 1D-output'
      endif
    endif

     if (index == first_index.and.parameters%id%id_type /= 'from_data') then
         open(unit, file = trim(fullname), status = 'replace', err = 120)
     else
         open(unit, file = trim(fullname), status = 'unknown', err = 120, &
    &         position = 'append')
     end if

    if (dim == 0) then  ! 0D - output
       write(unit, *) parameters%timer%t, yval
    else                ! 1D - output
       write(unit, *)
       write(unit, *)

       ! Write time information and then the actual (optionally downsampled)
       ! data. Add two empty lines to indicate the end of the time slice.
       write(tim,'(f15.10)') parameters%timer%t !for pygraph
       if (parameters%output%suffix.eq.'.yg') write(unit,'(a8,a)') '#Time = ', trim(adjustl(tim)) !'(a8,f16.10)' '(a8,e17.10e3)' '(a8,es12.5e3)'
!       if (parameters%output%suffix.eq.'.dat') write(unit,'(a8,f16.10)') '"Time = ', parameters%timer%t
       if (parameters%output%suffix.eq.'.dat') write(unit,'(a8,a)') '"Time = ', trim(adjustl(tim)) ! 2023/03/08: modified for compatibility with muninn
!write(*,*), first_point, last_point
!       do count = first_point + stagdiff(k) , last_point -coeflp/k, parameters%output%every_x_step
       do count = first_point + stagdiff(k) , last_point, parameters%output%every_x_step !change for staggered grid
          write(unit,'(2(e33.24e3))') xval(count), yval(count)
       end do
!       write(unit, *)
    endif
    close (unit) ! Close the output file.

    first = .false.
    return

120 call file_error(filename)
  end subroutine write_out_par

  subroutine file_error(filename)
    character(len = *), intent(in)    :: filename
    write (*, *) "File ", filename ," error: probably quota exceeded!"
    stop 'Program termination'
  end subroutine file_error


    subroutine write_t(k, index, xval, yval, filename, parameters)
    use parameters_mod
    implicit none
    integer,                           intent(in) :: index, k !!!!!!!!!!1
    real(kind=wp), 				       intent(in) :: xval, yval
    character(len=*),                  intent(in) :: filename
    type(parameters_type),             intent(in) :: parameters
    character(len=128)          :: fullname
    logical, save               :: first = .true.
    integer, parameter          :: unit = 1112
    integer, save               :: dim, first_index

    fullname = trim(file_prefix) // trim(filename) // trim(parameters%output%suffix)

    if (first) first_index = index
!print*, first

     if (index == first_index.and.parameters%id%id_type /= 'from_data') then
         open(unit, file = trim(fullname), status = 'replace')
     else
         open(unit, file = trim(fullname), status = 'unknown', position = 'append')
     end if

       write(unit, *) parameters%timer%t, yval

    close (unit) ! Close the output file.

    first = .false.
    return
  end subroutine write_t



  subroutine write_step(k, i, parameters, u, constr, x)!!!!!!!!!!!1
    use parameters_mod
    use numservicef90 ! added charac

    implicit none

    integer,                              intent(in)    :: i, k!!!!!!!!!!!!!!!111
    type(parameters_type),                intent(in)    :: parameters
    real(kind=wp), dimension(:,:),        intent(in)    :: u, constr
    real(kind=wp), dimension(:),          intent(in)    :: x

!    real(kind=wp)  :: out(size(u,1),size(x))
    real(kind=wp)  :: out(number_of_files,size(x))
    real(kind=wp)  :: dout(size(u,1),size(x)) ! added charac
    real(kind=wp)  :: charac(size(u,1),size(x)) ! added charac
    character :: num*2 ! added charac

    integer                     :: count

!changed

    out(1:size(u,1),:) = u
!    out(1:size(u,1),:) = log(abs(u))

!   out(1:size(u,1),:) = (u(:,2:size(u,2)+1)-u(:,0:size(u,2)-1))!/(x(2)-x(1))/2._wp ! first derivative
!   out(1:size(u,1),:) = (u(:,2:size(u,2)+1)-2*u(:,1:size(u,2))+u(:,0:size(u,2)-1))!/((x(2)-x(1))**2) ! second derivative
!   out(1:size(u,1),:) = (u(:,3:size(u,2)+2)-2*u(:,2:size(u,2)+1)+2*u(:,0:size(u,2)-1)-u(:,-1:size(u,2)-2))!/(2.0_wp*(x(2)-x(1))**3) ! third derivative
!out(1,:) = u(1,:)/x/x
!out(2,:) = u(2,:)/x/x
!out(3,:) = u(3,:)/x
!out(4,:) = u(4,:)/x/x
!out(5,:) = u(5,:)/x/x
!out(6,:) = x*u(6,:)
!out(7,:) = u(7,:)/x/x
!out(8,:) = u(8,:)/x
!out(9,:) = u(9,:)/x
!out(2,:) = (out(7,:)*out(5,:)+((u(7,2:size(u,2)+1)-u(7,0:size(u,2)-1))/(x(2)-x(1))/2._wp))/(2._wp)
!out(4,:) = (-out(7,:)*out(5,:)+((u(7,2:size(u,2)+1)-u(7,0:size(u,2)-1))/(x(2)-x(1))/2._wp))/(2._wp)

!out(7,:) = u(7,:)-(1+x*x)/2._wp
!out(8,:) = u(8,:)+x

!!! rescaling of phi -> phi*r/omega !!!
! out(10,:) = u(10,:)*2*parameters%physics%xscri*parameters%physics%aa*x/(parameters%physics%xscri-x)/(parameters%physics%xscri+x)
! if(parameters%grid%origin.eq.'nostag') out(10,size(x)-parameters%grid%nghost) = -parameters%physics%aa*parameters%physics%xscri*(u(10,size(x)-parameters%grid%nghost+1)-u(10,size(x)-parameters%grid%nghost-1))/(x(2)-x(1))/2._wp

!out(10,size(x)-parameters%grid%nghost)
!print*, x(size(x)-parameters%grid%nghost)

    out(size(u,1)+1:number_of_files,:) = constr

! added charac
  dout = u
  if (parameters%moldef%deriv_method == 'c8') then
     call mdiff_c8(dout(:,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
     call mdiff_c6(dout(:,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call mdiff_c4(dout(:,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
     call mdiff(dout(:,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if

!out(1,:) = dout(1,:)
!out(2,:) = dout(2,:)
!out(3,:) = dout(2,:)

  charac(1,:) = out(7,:)
  charac(2,:) = out(3,:)
  charac(3,:) = out(1,:)
  charac(4,:) = ( 2/3._wp*dout(7,:)*out(3,:)/out(7,:) - 2/3._wp*out(6,:)/out(3,:) - 1/3._wp*dout(1,:)*out(3,:)/out(1,:) )*x*2/(1-x**2)
  charac(5,:) = ( 2/3._wp*dout(7,:)*out(1,:)/out(7,:) + dout(1,:) )*x*2/(1-x**2)
  charac(6,:) = ( dout(7,:)/2._wp + out(5,:)*out(7,:)/(2._wp*out(3,:)*sqrt(out(1,:))) )*x*2/(1-x**2)
  charac(7,:) = ( dout(3,:)/2._wp - dout(7,:)*out(3,:)/3._wp/out(7,:) + out(6,:)/3._wp/out(3,:) + dout(1,:)*out(3,:)/6._wp/out(1,:) - out(4,:)*out(3,:)**2/2._wp/sqrt(out(1,:)) )*x*2/(1-x**2)
  charac(8,:) = ( dout(7,:)/2._wp - out(5,:)*out(7,:)/(2._wp*out(3,:)*sqrt(out(1,:))) )*x*2/(1-x**2)
  charac(9,:) = ( dout(3,:)/2._wp - dout(7,:)*out(3,:)/3._wp/out(7,:) + out(6,:)/3._wp/out(3,:) + dout(1,:)*out(3,:)/6._wp/out(1,:) + out(4,:)*out(3,:)**2/2._wp/sqrt(out(1,:)) )*x*2/(1-x**2)

       do count = 1, 9
    write(num,'(i2)') count
    num = adjustl(num)
!    call write_out_par (k, i, x, charac(count,:), "char"//trim(num)//'_'//trim(parameters%moldef%deriv_method), parameters)
       end do
! added charac

!        do count = 1, size(u, 1)
        do count = 1, number_of_files
    ! I am outputting only phi now. If I want the others, I shall replace u(3,:) with u(count,:) and file_names(1) with file_names(count) (2 May 2008)
    call write_out_par (k, i, x, out(count,:), file_names(count), parameters)!!!!!!!!!!!!!!!!!!!!!!!!!11
        end do
!changed


  end subroutine write_step

    subroutine change_file_names(file_names, parameters, string)

      type(parameters_type),                            intent(in)    :: parameters
      character(len = 120), dimension(number_of_files), intent(inout) :: file_names
      character(len = *),                               intent(in)    :: string
      integer i

      do i = 1, number_of_files
         file_names(i) = trim(file_names(i)) // '_' // &
              trim(parameters%moldef%deriv_method) // '_' // &
              string // ''
    enddo
  end subroutine change_file_names

end module output
