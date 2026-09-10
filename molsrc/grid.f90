
module grid
  use prec
  use parameters_mod
  use molparameters

  implicit none

  real(kind = iwp), dimension(:), allocatable :: gridvect

contains

subroutine set_default_grid_parameters
    implicit none
    default_parameters%grid%ncells    =  99
    default_parameters%grid%npoints   = 100
    default_parameters%grid%nghost    =   0
    default_parameters%grid%xmin      =  -0.5_wp
    default_parameters%grid%xmax      =   0.5_wp
    default_parameters%grid%dx        =   0.01_wp
    default_parameters%grid%fill_type =   'undefd'
    default_parameters%grid%origin    =   'nostag'
  end subroutine set_default_grid_parameters

  subroutine grid_param_testing(parameters)
    implicit none

    type(parameters_type), intent(inout) :: parameters

    parameters%grid%npoints = parameters%grid%ncells + 1
    if(parameters%grid%ncells < 0) then
       stop 'ncells too small!'
    endif

    if(parameters%grid%nghost < 0) then
       stop 'nghost too small!'
    endif

    if(parameters%grid%xmax < parameters%grid%xmin .and. &
         parameters%grid%fill_type == 'minmax') then
       stop 'xmin and/or xmax are wrong!'
    endif

    if(parameters%grid%dx <=0 .and. parameters%grid%fill_type == 'mindx') then
       stop 'dx too small!'
    endif
   end subroutine grid_param_testing


  subroutine grid_allocate(lbnd, ubnd)
    integer, intent(in) :: lbnd, ubnd

    allocate(gridvect(lbnd : ubnd))
  end subroutine grid_allocate

  subroutine grid_deallocate
    deallocate(gridvect)  
  end subroutine grid_deallocate


  subroutine new_grid(parameters)

    type(parameters_type),       intent(in) :: parameters

    integer                     :: npoints, nghost, allpoints, ncells, counter
    real(kind = wp)             :: dx, xmax, xmin

    ! executable statements

    allpoints = size(gridvect)
    nghost    = lbound(gridvect, 1)

    if (nghost <= 0) then
       nghost = abs(nghost) 
    else
       write(*,*) nghost
       STOP  'inconsistent lower bound of grid [has to be <= 0]'
    endif

    npoints   = allpoints - 2 * nghost 
    ncells    = npoints   - 1 


    do counter = lbound(gridvect, 1), ubound(gridvect, 1)
       gridvect(counter) = counter
    enddo


    xmin = parameters%grid%xmin
    if(ncells == 0 .and. parameters%grid%fill_type == "minmax") STOP 'ncells=0 and fill_type=minmax are not compatible'
    if (parameters%grid%fill_type == "minmax" .and. ncells .ne. 0) then
       xmax = parameters%grid%xmax
       dx   = abs(xmax - xmin) / dble(ncells)
       if (parameters%grid%origin.eq.'stag') dx   = abs(xmax - xmin) / (dble(ncells+1)) !change for staggered grid
       if (parameters%grid%origin.eq.'misman') dx   = abs(xmax - xmin) / (dble(ncells)+0.5_wp) !06-07-2016 misman
    elseif (parameters%grid%fill_type == "min_dx") then
       dx = parameters%grid%dx
    else
       STOP  'filltype has to be either "minmax" or "min_dx"'
    endif

    if (parameters%grid%origin.eq.'nostag') then
	gridvect = xmin + gridvect * dx                
    elseif (parameters%grid%origin.eq.'stag') then
	gridvect = xmin + gridvect * dx                    + dx/2.0_wp  !!!!!!!staggered grid                
    elseif (parameters%grid%origin.eq.'misman') then
	gridvect = xmin + gridvect * dx                    + dx/2.0_wp  !06-07-2016 misman
    end if

    write(*, *) 'grid spacing           =  ', dx
    write(*, *) 'grid min, max          =  ', gridvect(0), gridvect(ncells)
    write(*, *) 'ghost zones            =  ', nghost
    if(dx>0)    write(*, *) 'relative fill mismatch =  ', abs(gridvect(ncells) - xmax)/dx

!    open(10, file='grid.xg')
!    do counter = lbound(gridvect, 1), ubound(gridvect, 1)
!       write (10,*) counter, gridvect(counter)
!    enddo
!    close(10)



  end subroutine new_grid
end module grid

