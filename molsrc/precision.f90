module prec
  implicit none
  ! Precision parameter
!  integer, parameter :: wp = selected_real_kind(p=32)
  integer, parameter :: wp = 8 !kind(1.0d0)
  integer, parameter :: iwp = 16 !kind(1.0d0)
  real(kind=wp), parameter :: pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068d0
end module prec
