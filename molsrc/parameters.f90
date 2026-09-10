module parameters_mod

 use molparameters
 use userparameters

 implicit none

  ! The definitions of the different parameter types
  type parameters_type
     type(physics_type) :: physics
     type(id_type)      :: id
     type(slice_type)   :: slice
     type(time_type)    :: timer
     type(moldef_type)  :: moldef
     type(grid_type)    :: grid
     type(output_type)  :: output
  end type parameters_type

  type(parameters_type) :: default_parameters

end module parameters_mod
