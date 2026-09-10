module molparameters

  use prec
  implicit none

 type moldef_type
     character(len = 42):: boundary_type
     character(len = 42):: stencil_boundary !to consider one-sided stencils at the boundaries
     character(len = 42):: eqsystem !linear or nonlinear system
     real(kind = wp)    :: dissipation_eps  !value of epsilon of the dissipation
     character(len = 42):: spatial_order
     character(len = 42):: int_method
     character(len = 42):: deriv_method
  end type moldef_type

  type time_type
     real(kind = wp)    :: dt          !time step to the new slice
     real(kind = wp)    :: t           !current time
     real(kind = wp)    :: end_time    !time to stop evolution
     integer            :: step        !current time step
  end type time_type

  type grid_type
     integer            :: ncells
     integer            :: npoints
     integer            :: nghost
     real(kind = wp)    :: xmin
     real(kind = wp)    :: xmax
     real(kind = wp)    :: dx
     real(kind = wp)    :: boundary_factor
     character(len = 6) :: fill_type
     character(len = 6) :: origin !option to stagger the grid at the origin or not
     character(len = 6) :: originbc ! boundary condition at origin: parity for regular spacetime or extrapolation for trumpet or excision
  end type grid_type

 type output_type
     integer            :: every_t_step
     integer            :: every_x_step
     character(len = 42):: suffix
  end type output_type

  ! Input - Output parameters

  integer, parameter :: number_of_files = 28
  character(len = 120), dimension(number_of_files) :: file_names, file_names_old
  character(len = 120) :: file_prefix


end module molparameters
