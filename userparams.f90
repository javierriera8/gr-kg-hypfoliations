module userparameters

  use prec

  implicit none

  ! The definitions of the different parameter types
  type physics_type
     character(len = 42):: equation
     real(kind=wp) 	:: v, eta, lambda, mu, mass, cosm, aa, k1, k2, xscri, coealphah, coealphal, coebetar, coeLambdar, coeZr, pphi, trafK, Kcmc, Ccmc, propsign, Q, xexc, nuplogv, nconstv, masskgf, qp
     integer 		:: kgf, z4, z4c, omt, em, gr
  end type physics_type

  ! For the initial data Gauss function
  type id_type
     character(len = 42):: id_type, simulation ! simulation can be sequencing, bisection or single
     real(kind=wp) :: center, a, sigma, da, centero, ao, sigmao, aim, sigmaim, centerim ! da is for sequencing runs.
     integer		:: noise
  end type id_type

  ! For the slicing
  type slice_type
     real(kind=wp) :: cuplog, charm, rmatchmin, rmatchmax, h, tvaldiag, tcol !, k
     integer	   :: evollapse, evolshift
     ! h... cmc_ss, the integration constant of cmc surface
     ! k... cmc_ss, mean curvature
  end type slice_type


!  real(kind=wp), parameter :: pi = 3.141592653589793238462643383279

end module userparameters
