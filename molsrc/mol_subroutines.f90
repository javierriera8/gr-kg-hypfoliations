!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines that were originally located in the module mol,
! moved here for convenience. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mol_subroutines
  use prec
  use parameters_mod
  use molparameters

  implicit none

contains


!!!!!!!!!!!!!!!!!! Setting the expression for the potential of the scalar field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! derivative of the potential with respect to the scalar field
subroutine Vkgprimeset(x,parameters,omega,rephi,imphi,Vkgprimefunc)
implicit none
type(parameters_type), intent(inout)  :: parameters
real(kind=wp), dimension(:),   intent(in)  :: x
real(kind=wp), dimension(size(x)) :: omega, rephi, imphi, Vkgprimefunc
real(kind=wp) :: masskgf

masskgf = parameters%physics%masskgf

Vkgprimefunc = (1d0/2)*(masskgf**2) - (1d0/2)*(100)*(omega**2)*(rephi**2 + imphi**2) + (1d0/2)*(200)*(omega**4)*(rephi**2 + imphi**2)**2

!POTENTIALS PRIME TESTED
!potential1 = (masskgf**2)
!potential2 = (masskgf**2)*(omega**2)*2*(rephi**2 + imphi**2)
!potential3 = (masskgf**2)*(omega**4)*3*(rephi**2 + imphi**2)**2
!potential4 = (1d0/2)*(masskgf**2) - (1d0/2)*(10240d0)*(omega**2)*(rephi**2 + imphi**2) + (1d0/2)*(2.18d0 * 10**7)*(omega**4)*(rephi**2 + imphi**2)**2
!potential5 = (masskgf**2)*cos((omega**2)*(rephi**2 + imphi**2))

return
end subroutine Vkgprimeset

! potential
subroutine Vkgset(x,parameters,omega,rephi,imphi,Vkgfunc)
implicit none
type(parameters_type), intent(inout)  :: parameters
real(kind=wp), dimension(:),   intent(in)  :: x
real(kind=wp), dimension(size(x)) :: omega, rephi, imphi, Vkgfunc
real(kind=wp) :: masskgf!, rmatchmax
!integer :: imatchmax

masskgf = parameters%physics%masskgf
!rmatchmax = parameters%slice%rmatchmax
                                                                          !lambda 10240d0                                            !mu 2.18d7 
Vkgfunc =(1d0/2)*(masskgf**2)*(omega**2)*(rephi**2 + imphi**2) - (1d0/4)*(100)*(omega**4)*(rephi**2 + imphi**2)**2 + (1d0/6)*(200)*(omega**6)*(rephi**2 + imphi**2)**3

!change the potential at: Vkgprimeset, Vkgset and Vkgfval
!POTENTIALS TESTED
!potential1 = m^2 phi^2        (masskgf**2)*(omega**2)*(rephi**2 + imphi**2)
!potential2 = m^2 (phi^2)^2    (masskgf**2)*(omega**4)*(rephi**2 + imphi**2)**2
!potential3 = m^2 (phi^2)^3    (masskgf**2)*(omega**6)*(rephi**2 + imphi**2)**3
!potential4 = 0.5m^2 phi^2 - 0.25 lambda (phi^2)^2 + (1/6) mu (phi^2)^3 !we test lambda=10240 i mu=2.18*107    (1d0/2)*(masskgf**2)*(omega**2)*(rephi**2 + imphi**2) - (1d0/4)*(10240d0)*(omega**4)*(rephi**2 + imphi**2)**2 + (1d0/6)*(2.18d0 * 10**7)*(omega**6)*(rephi**2 + imphi**2)**3
!potential5 = m^2 sin(phi^2)   (masskgf**2)*sin((omega**2)*(rephi**2 + imphi**2))

!* (1 - x**2)**4 !function8 ! change here and in following function masskgfval

!below function9
!masskgfunc = 0._wp
!imatchmax = minloc(abs(x-rmatchmax),1) ! find closest point to outer matching point
!masskgfunc(:imatchmax) = masskgf*( 1 - (70*x(:imatchmax)**9)/rmatchmax**9 + (315*x(:imatchmax)**8)/rmatchmax**8 - (540*x(:imatchmax)**7)/rmatchmax**7 + (420*x(:imatchmax)**6)/rmatchmax**6 - (126*x(:imatchmax)**5)/rmatchmax**5 ) ! simple suggestion obtained by chatgpt when asking for up to 4 continuous derivatives

! it may be worth to try with a piecewise function that is exactly flat in the central part and then decreases towards scri

!FUNCTIONS TESTED
!function1 = 1 - x**2
!function2 = (1 - x**2)**2
!function3 = cos((pi * x)/2) !we have to use pi= 3.14159265358979323846
!function4 = exp(- (x**2/(1 - x**2))) !we have to define it as 0 where x >= 1 to avoid diving by 0
!function5 =  0.5d0 * ( 1.0d0 + tanh( (x - 0.8d0) / 0.05d0 ) ) funció esglaó suavitzada
!function6 = (1 - x**4)**2
!function7 = (1 - x**2)**3
!function8 = (1 - x**2)**4
return
end subroutine Vkgset

! potential -- used in the intitial data calculations -- need to generalize to include non-zero imphi in inidata as well
function Vkgfval(x, masskgf, omega, rephiphys, imphiphys)
implicit none
real(kind=iwp) :: x, masskgf, Vkgfval, omega, rephiphys, imphiphys
  Vkgfval = (1d0/2)*(masskgf**2)*(rephiphys**2 + imphiphys**2) - (1d0/4)*(100)*(rephiphys**2 + imphiphys**2)**2 + (1d0/6)*(200)*(rephiphys**2 + imphiphys**2)**3
  !no omega, as here we are including the physical scalar field, for the way the initial data solver for the Hamiltonian constraint is set
  
   !* (1 - x**2)**4 !function8 ! change here and in previous subroutine Vkgfset

 !POTENTIALS TESTED
 !potential1 = (masskgf**2)*(rephiphys**2 + imphiphys**2)
 !potential2 = (masskgf**2)*(rephiphys**2 + imphiphys**2)**2 
 !potential3 = (masskgf**2)*(rephiphys**2 + imphiphys**2)**3 
 !potential4 = (1d0/2)*(masskgf**2)*(rephiphys**2 + imphiphys**2) - (1d0/4)*(10240d0)*(rephiphys**2 + imphiphys**2)**2 + (1d0/6)*(2.18d0 * 10**7)*(rephiphys**2 + imphiphys**2)**3
 !potential5 = (masskgf**2)*sin(rephiphys**2 + imphiphys**2)



  !below function9
!  masskgfval = 0._iwp
 ! if (x.lt.rmatchmax) masskgfval = masskgf*( 1 - (70*x**9)/rmatchmax**9 + (315*x**8)/rmatchmax**8 - (540*x**7)/rmatchmax**7 + (420*x**6)/rmatchmax**6 - (126*x**5)/rmatchmax**5 )
  return
end function Vkgfval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tests for dynamical calculation of the gauge source terms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!! Ccmc calculation !!!!!!!!!
function ccmccal(Kcmc,mass,Q)
implicit none
real(kind=iwp), intent(in)  :: Kcmc, mass, Q
real(kind=iwp) :: icond1, icond2, icond3, epsi1, epsi2, epsi3
real(kind=iwp) :: ccmccal
integer :: nt
    icond1 = 0
    icond2 = 20*(-Kcmc)*mass + 1.5_iwp*mass ! experimental - may have to be changed
!    if (abs(Q-mass).lt.5d-2) icond2 = -Kcmc*mass**3/3._iwp + 50*abs(mass-Q)**0.7_iwp*mass ! problems with Q almost mass, as curve beyond Ccmc solution is very steep and initial condition is too large (negatively)
  nt = 0
  epsi1 = ccmcval(Kcmc,mass,Q,icond1)
  epsi2 = ccmcval(Kcmc,mass,Q,icond2)
  648 continue
  nt = nt + 1
!  print*, nt, epsi3
  if (epsi1*epsi2.gt.0) then
  write(*,*) "Initial values for the bisection for Ccmc not correct: "
  write(*,'(4(f10.4,a))') icond1, ' -> ', epsi1, ',  ', icond2, ' -> ', epsi2
  stop
  endif
  icond3 = (icond2+icond1)/2._iwp
  epsi3 = ccmcval(Kcmc,mass,Q,icond3)
  if (abs(epsi3).gt.1d-22) then !12
    if (epsi1*epsi3.gt.0) then
 icond1 = icond3
 epsi1 = epsi3
    else
 icond2 = icond3
 epsi2 = epsi3
    endif
    if (nt.gt.100000) stop "Arrived at 100000 iterations for Ccmc - decrease accuracy."
    goto 648
  endif
  ccmccal = icond3
return
end function ccmccal

function ccmcval(Kcmc,mass,Q,Ccmc)
implicit none
real(kind=iwp), intent(in)  :: Kcmc, mass, Q, Ccmc
real(kind=iwp) :: ccmcval
ccmcval = (-64*(1296*Ccmc**4 + 729*Ccmc**6*Kcmc**4 - 2916*Ccmc**5*Kcmc**3*mass + 3402*Ccmc**4*Kcmc**2*mass**2 + 2916*Ccmc**3*Kcmc*mass**3 + 1944*Ccmc**5*Kcmc**5*mass**3 - 2187*Ccmc**2*mass**4 - 8748*Ccmc**4*Kcmc**4*mass**4 + 13122*Ccmc**3*Kcmc**3*mass**5 - 6561*Ccmc**2*Kcmc**2*mass**6 - 972*Ccmc**4*Kcmc**2*Q**2 - 1944*Ccmc**3*Kcmc*mass*Q**2 - 1458*Ccmc**5*Kcmc**5*mass*Q**2 + 2916*Ccmc**2*mass**2*Q**2 + 7047*Ccmc**4*Kcmc**4*mass**2*Q**2 - 14580*Ccmc**3*Kcmc**3*mass**3*Q**2 + 10935*Ccmc**2*Kcmc**2*mass**4*Q**2 - 648*Ccmc**2*Q**4 - 243*Ccmc**4*Kcmc**4*Q**4 + 3078*Ccmc**3*Kcmc**3*mass*Q**4 - 4617*Ccmc**2*Kcmc**2*mass**2*Q**4 + 603*Ccmc**2*Kcmc**2*Q**6 + 81*Ccmc**4*Kcmc**6*Q**6 + 54*Ccmc*Kcmc*mass*Q**6 - 324*Ccmc**3*Kcmc**5*mass*Q**6 - 81*mass**2*Q**6 + 378*Ccmc**2*Kcmc**4*mass**2*Q**6 + 324*Ccmc*Kcmc**3*mass**3*Q**6 - 243*Kcmc**2*mass**4*Q**6 + 81*Q**8 - 108*Ccmc**2*Kcmc**4*Q**8 - 216*Ccmc*Kcmc**3*mass*Q**8 + 324*Kcmc**2*mass**2*Q**8 - 72*Kcmc**2*Q**10 + 16*Kcmc**4*Q**12))/729._iwp !Ccmc**2*Kcmc**2* solution goes to zero if included ! still does not work for extreme Q=mass cas, as solution is a maximum -> using first derivative for this case
if (abs(Q-mass).lt.1d-12) ccmcval = (-128*Kcmc**2*(243*Ccmc**5*Kcmc**4 + 3*Kcmc*mass*Q**6*(1 + Kcmc**2*(6*mass**2 - 4*Q**2)) + 135*Ccmc**4*Kcmc**3*mass*(-6 + Kcmc**2*(4*mass**2 - 3*Q**2)) + 18*Ccmc**3*(16 + Kcmc**6*Q**6 + 6*Kcmc**2*(7*mass**2 - 2*Q**2) - 3*Kcmc**4*(36*mass**4 - 29*mass**2*Q**2 + Q**4)) - 27*Ccmc**2*Kcmc*mass*(-81*Kcmc**2*mass**4 + 18*mass**2*(-1 + 5*Kcmc**2*Q**2) + Q**2*(12 - 19*Kcmc**2*Q**2 + 2*Kcmc**4*Q**4)) + Ccmc*(-729*Kcmc**2*mass**6 - 72*Q**4 + 67*Kcmc**2*Q**6 - 12*Kcmc**4*Q**8 + 243*mass**4*(-1 + 5*Kcmc**2*Q**2) + mass**2*(324*Q**2 - 513*Kcmc**2*Q**4 + 42*Kcmc**4*Q**6))))/81._iwp
if (abs(Q).lt.1d-12) ccmcval = (-64*(9*Ccmc**4*Kcmc**4 + 12*Ccmc**3*Kcmc**3*mass*(-3 + 2*Kcmc**2*mass**2) + 18*Ccmc*Kcmc*mass**3*(2 + 9*Kcmc**2*mass**2) + Ccmc**2*(16 + 42*Kcmc**2*mass**2 - 108*Kcmc**4*mass**4) - 27*(mass**4 + 3*Kcmc**2*mass**6)))/9._iwp !Ccmc**4*Kcmc**2* solution goes to zero if included ! if Q=0, then discriminant has a Ccmc**4
return
end function ccmcval

!!!!!!!!! inmstR calculation !!!!!!!!!

function inmstRcal(Kcmc,mass,Q,Ccmc)
implicit none
real(kind=iwp), intent(in)  :: Kcmc, mass, Q, Ccmc
real(kind=iwp) :: icond1, icond2, icond3, epsi1, epsi2, epsi3
real(kind=iwp) :: inmstRcal
integer :: nt
    icond1 = 0.9_iwp*mass
    icond2 = 2*mass
  nt = 0
  epsi1 = inmstRval(Kcmc,mass,Q,Ccmc,icond1)
  epsi2 = inmstRval(Kcmc,mass,Q,Ccmc,icond2)
  649 continue
  nt = nt + 1
!  print*, nt, epsi3
  if (epsi1*epsi2.gt.0) then
  write(*,*) "Initial values for the bisection for inmstR not correct: "
  write(*,'(4(f10.4,a))') icond1, ' -> ', epsi1, ',  ', icond2, ' -> ', epsi2
  stop
  endif
  icond3 = (icond2+icond1)/2._iwp
  epsi3 = inmstRval(Kcmc,mass,Q,Ccmc,icond3)
  if (abs(epsi3).gt.1d-22) then !12
    if (epsi1*epsi3.gt.0) then
 icond1 = icond3
 epsi1 = epsi3
    else
 icond2 = icond3
 epsi2 = epsi3
    endif
    if (nt.gt.100000) stop "Arrived at 100000 iterations for inmstR - decrease accuracy."
    goto 649
  endif
  inmstRcal = icond3
return
end function inmstRcal

function inmstRval(Kcmc,mass,Q,Ccmc,x)
implicit none
real(kind=iwp), intent(in)  :: Kcmc, mass, Q, Ccmc, x
real(kind=iwp) :: inmstRval
inmstRval = 18*Q**2*x + 18*(Ccmc*Kcmc - 3*mass)*x**2 + 36*x**3 + 6*Kcmc**2*x**5 ! first derivative
!(18*Ccmc*Kcmc - 54*mass)*x**2 + 36*x**3 + 6*Kcmc**2*x**5 ! first derivative
!9*Ccmc**2 + (6*Ccmc*Kcmc - 18*mass)*x**3 + 9*x**4 + Kcmc**2*x**6 !function, not appropriate, as the double root is a minimum and cannot be calculated using the bisection method (plot the polynomial to understand). however, first derivative has the same root (simple in this case).
return
end function inmstRval

!!!!!!! Calculation of new psic!!!!!!!!!!

subroutine dynpsic(x,parameters,psic)
implicit none
integer, parameter :: fac = 1
type(parameters_type), intent(inout)  :: parameters
real(kind=iwp), dimension(:),   intent(inout)  :: x
real(kind=iwp) :: mass, Kcmc, Q, Ccmc, aa, xscri, dx
real(kind=wp), dimension(size(x)) :: psic
real(kind=wp), dimension(size(x)) :: omega
real(kind=iwp), dimension(size(x)*fac*2) :: psicompl !, allocatable
integer :: counter

mass = parameters%physics%mass
Kcmc = parameters%physics%Kcmc
Ccmc = parameters%physics%Ccmc
Q = parameters%physics%Q
aa = parameters%physics%aa
xscri = parameters%physics%xscri

if (aa.lt.0) aa = -3/Kcmc

dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / real(parameters%grid%ncells,iwp)

!print*, Ccmc, Kcmc, mass, Q
Ccmc = ccmccal(Kcmc,mass,Q)
parameters%physics%Ccmc = Ccmc !! dangerous??
print*, "New Ccmc = ", Ccmc
!print*, mass, Ccmc, Kcmc

omega = (xscri+x)*(xscri-x)/(2._wp*aa*xscri)

if((abs(mass).gt.1d-12).and.(abs(Ccmc).gt.1d-12))then
  write(*,*) "Determination of the overall compactification factor:"
  call aconfcal(fac,x,dx,parameters,psic,psicompl)
	do counter = 1, abs(parameters%grid%nghost) ! bcs for aconf(=psic)
		psic(lbound(x, 1) + counter - 1) = - psic(lbound(x, 1) + 2*abs(parameters%grid%nghost) - counter)
		psic(ubound(x, 1)-counter+1) = omega(ubound(x, 1)-counter+1)
	enddo
  do counter = lbound(x, 1), ubound(x, 1)
	write(566,*) x(counter), psic(counter), omega(counter)
  enddo
else
write(*,*) "As mass = 0 and Ccmc = 0, aconf will be set to omega."
  psic = omega
endif

return
end subroutine dynpsic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculation of the overall compactification factor (called here aconf and psic in the rest of the code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aconfcal(fac,x,dx,parameters,psic,aconfc)
  implicit none
  real(kind=iwp), dimension(:),   intent(in)  :: x
  type(parameters_type),intent(in)  :: parameters
  real(kind=wp), dimension(size(x)) :: psic
  real(kind=iwp) :: eaconf1, eaconf2, eaconf3, icond1, icond2, icond3
  real(kind=iwp) :: dx
  integer, parameter :: facnostag=2 ! do not change, to calculate inital data with nostag grid
  integer :: lb, li, ui, ub, nt, fac !, counter
  real(kind=iwp), dimension(:) :: aconfc !, allocatable

!    allocate(aconfc(1:size(x)*fac*facnostag))

    !Original number of points
    lb = lbound(x,1)
    ub = ubound(x,1)
    li = lb + parameters%grid%nghost
    ui = ub - parameters%grid%nghost
!print*, lb, li, ui, ub, dx
    !Increased number of points (to calculate with higher resolution)
    dx=dx/fac/facnostag
    ui = li + parameters%grid%ncells*fac*facnostag - 1
    if(parameters%grid%origin.eq.'nostag'.or.parameters%grid%origin.eq.'misman') ui = li + parameters%grid%ncells*fac*facnostag
    ub = ui + parameters%grid%nghost
!print*, lb, li, ui, ub, dx
!-!   if (parameters%grid%origin.eq.'nostag') stop "Calculation of aconf not prepared for a non-staggered grid."  !-!

  icond1 = 2000._iwp !600._iwp ! 300._iwp !0.05_iwp !1000K=-0.5 !1.5_wpK=-1 !0.5_wpK=-3 !0.3_wp !0.9_wp
  icond2 = -0.0002_iwp !0.7_wp !1.0_wp
  nt = 0
  if(abs(parameters%grid%xmin).gt.1d-12) stop "Leftmost boundary not set at the imposed xmin = 0."
  call eaconfc(icond1,lb,li,ui,ub,0._iwp+dx,dx,eaconf1,aconfc,parameters)
!print*, icond1, eaconf1
  call eaconfc(icond2,lb,li,ui,ub,0._iwp+dx,dx,eaconf2,aconfc,parameters)
!print*, icond2, eaconf2
  659 continue
  nt = nt + 1
!  print*, eaconf1, eaconf2
! write(555,*) icond1, eaconf1
! write(555,*) icond2, eaconf2
  if (eaconf1*eaconf2.gt.0) then
  write(*,*) "Initial values for the bisection of aconf not correct: "
  write(*,'(4(f10.4,a))') icond1, ' -> ', eaconf1, ',  ', icond2, ' -> ', eaconf2
  stop
  endif
  icond3 = (icond2+icond1)/2._wp
  call eaconfc(icond3,lb,li,ui,ub,0._iwp+dx,dx,eaconf3,aconfc,parameters)
!print*, icond3, eaconf3
!stop
  if (abs(eaconf3).gt.1d-20) then !12
    if (eaconf1*eaconf3.gt.0) then
 icond1 = icond3
 eaconf1 = eaconf3
    else
 icond2 = icond3
 eaconf2 = eaconf3
    endif
    goto 659
  endif
  write(*,*) '-- number of iterations = ', nt, eaconf3
  write(*,*) 'selected initial condition = ', icond3

    !Back to original number of points
    dx=dx*fac*facnostag
    ui = li + parameters%grid%ncells - 1
    if(parameters%grid%origin.eq.'nostag'.or.parameters%grid%origin.eq.'misman') ui = li + parameters%grid%ncells
    ub = ui + parameters%grid%nghost
!print*, lb, li, ui, ub, dx
    psic(li:ui) = aconfc(li+int((fac-1)/2._wp)+int(facnostag/2._wp):ui:fac*facnostag) ! stag and misman
	if(parameters%grid%origin.eq.'nostag') psic(li:ui) = aconfc(li:ui:fac*facnostag) !-! ! nostag
!	do counter = lbound(aconfc,1), ubound(aconfc,1)
!		print*, counter, aconfc(counter), li, ui
!	enddo
!stop
  return
  end subroutine aconfcal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eaconfc(icond,lb,li,ui,ub,x,dx,eaconf,aconfc,parameters)
  implicit none
  type(parameters_type), intent(in)  :: parameters
  real(kind=iwp), intent(in)  :: x
  real(kind=iwp) :: eaconf, icond
  real(kind=iwp) :: mass, Kcmc, Q, Ccmc, inmstR
  real(kind=iwp) :: dx, xv
  integer :: lb, li, ui, ub, counter, maxcount
  real(kind=iwp) :: aconf, aconfm, arhs, k1, k2, k3, k4
  real(kind=iwp), dimension(lb:ub) :: aconfc

  mass = parameters%physics%mass
  Kcmc = parameters%physics%Kcmc
  Ccmc = parameters%physics%Ccmc
  Q = parameters%physics%Q

inmstR = inmstRcal(Kcmc,mass,Q,Ccmc)
!print*, Kcmc, Ccmc, inmstR
  xv = x

!  aconf = icond+xv*(icond/xv - sqrt(9*Ccmc**2*icond**6+6*(Ccmc*Kcmc-3*mass)*icond**3*xv**3+9*icond**2*xv**4+Kcmc**2*xv**6)/(3._wp*xv**3))
!  aconf = icond+ icond*(inmstR**4+inmstR**3*(Ccmc*Kcmc-3*mass)*icond+3*Ccmc**2*icond**4)/inmstR**4 ! Taylor order 2
!  aconf = icond+ (-sqrt(inmstR**6*Kcmc**2+9*inmstR**4*icond**2+6*inmstR**3*(Ccmc*Kcmc-3*mass)*icond**3+9*Ccmc**2*icond**6)/(3._wp*inmstR**2)) ! Taylor order 1
  aconf = dx**2*icond
!print*, xv, aconf
write(555,*) xv, aconf !*xv/inmstR
  if(parameters%grid%origin.eq.'nostag') aconfc(li) = 0._iwp ! set aconf to zero at the origin
  aconfc(li+1) = -(aconf-1)*xv/inmstR
  maxcount = ui+1 ! stag and misman
  if(parameters%grid%origin.eq.'nostag') maxcount = ui ! nostag
  do counter = li+2, maxcount !ub!ui ! for staggered grid I am using ui+1
  aconfm = aconf
  call a2rhs(xv,aconfm,mass,Kcmc,Q,Ccmc,inmstR,arhs)
  if(.true.)then
!RK 4th order
  k1 = dx*arhs
  xv = xv+dx/2._iwp
  call a2rhs(xv,aconfm+k1/2._iwp,mass,Kcmc,Q,Ccmc,inmstR,arhs)
  k2 = dx*arhs
  call a2rhs(xv,aconfm+k2/2._iwp,mass,Kcmc,Q,Ccmc,inmstR,arhs)
  k3 = dx*arhs
  xv = xv+dx/2._iwp
  call a2rhs(xv,aconfm+k3,mass,Kcmc,Q,Ccmc,inmstR,arhs)
  k4 = dx*arhs
  aconf = aconfm + (k1+2*k2+2*k3+k4)/6._iwp
  else
!1st order
  aconf = aconfm + dx*arhs
  xv = xv + dx
  endif
  aconfc(counter) = -(aconf-1)*xv/inmstR
!print*, xv, aconf
write(555,*) xv, aconf !*xv/inmstR
  if (abs(aconf-aconfm).gt.1) aconf = aconfm ! keep??
  end do
!erpsi = psi(1)-1 ! attempt for non-staggered
  eaconf = aconf -1._iwp !(aconf+aconfm)/2._wp !psi(1)**4 - 1+4*pi/3._wp*dphikg*dphikg !psi(2) ! condition ! fix
!write(555,*)
  return
  end subroutine eaconfc

  subroutine a2rhs(x,aconfm,mass,Kcmc,Q,Ccmc,inmstR,arhs)
  implicit none
  real(kind=iwp) :: x, mass, Kcmc, Q, Ccmc, inmstR
  real(kind=iwp) :: aconfm, arhs
! no  arhs = aconfm/x - sqrt(9*Ccmc**2*aconfm**6+6*(Ccmc*Kcmc-3*mass)*aconfm**3*x**3+9*aconfm**2*x**4+Kcmc**2*x**6)/(3._iwp*x**3)
!  arhs = -sqrt(inmstR**6*Kcmc**2+9*inmstR**4*aconfm**2+6*inmstR**3*(Ccmc*Kcmc-3*mass)*aconfm**3+9*Ccmc**2*aconfm**6)/(3._iwp*inmstR**2*x) ! eaconfc
  arhs = Sqrt(9*Ccmc**2 + 6*Ccmc*inmstR**3*Kcmc + inmstR**2*(9*inmstR**2 + inmstR**4*Kcmc**2 - 18*inmstR*mass + 9*Q**2) - 18*(3*Ccmc**2 + Ccmc*inmstR**3*Kcmc + inmstR**2*(inmstR**2 - 3*inmstR*mass + 2*Q**2))*aconfm + 9*(15*Ccmc**2 + 2*Ccmc*inmstR**3*Kcmc + inmstR**2*(inmstR**2 - 6*inmstR*mass + 6*Q**2))*aconfm**2 - 6*(30*Ccmc**2 + Ccmc*inmstR**3*Kcmc - 3*inmstR**3*mass + 6*inmstR**2*Q**2)*aconfm**3 + 9*(15*Ccmc**2 + inmstR**2*Q**2)*aconfm**4 - 54*Ccmc**2*aconfm**5 + 9*Ccmc**2*aconfm**6)/(3._iwp*inmstR**2*x)
!  if (abs(x).lt.10d-12) arhs = 0._iwp !-!
!  arhs = Sqrt(9*Ccmc**2 + 135*aconfm**4*Ccmc**2 - 54*aconfm**5*Ccmc**2 + 9*aconfm**6*Ccmc**2 + 6*Ccmc*inmstR**3*Kcmc + 9*aconfm**2*(15*Ccmc**2 + 2*Ccmc*inmstR**3*Kcmc + inmstR**3*(inmstR - 6*mass)) - 18*aconfm*(3*Ccmc**2 + Ccmc*inmstR**3*Kcmc + inmstR**3*(inmstR - 3*mass)) + inmstR**3*(9*inmstR + inmstR**3*Kcmc**2 - 18*mass) - 6*aconfm**3*(30*Ccmc**2 + Ccmc*inmstR**3*Kcmc - 3*inmstR**3*mass))/(3._iwp*inmstR**2*x) ! for aconf->-(aconf-1) ! working one for Schwarzschild
!  arhs = aconfm/x - aconfm*sqrt(Ccmc**2*aconfm**4-2*mass*aconfm*x**3+x**4)/(x**3) ! eaconfcback
  return
  end subroutine a2rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!! calculation of rbaralpha, the innermost radius in the pseudo-compactified domain achievable by solving a choice for the cK slicing condition !!!!!!!!!
function rbaralphacal(Kcmc,mass,Q,xscri,coealphal,coealphah)
implicit none
real(kind=iwp), intent(in)  :: Kcmc, mass, Q, xscri, coealphal, coealphah
real(kind=iwp) :: icond1, icond2, icond3, epsi1, epsi2, epsi3
real(kind=iwp) :: rbaralphacal
integer :: nt
    icond1 = 0.1_iwp
    icond2 = 0.9_iwp*xscri
  nt = 0
  epsi1 = rbaralphaval(Kcmc,mass,Q,xscri,coealphal,coealphah,icond1)
  644 continue
  epsi2 = rbaralphaval(Kcmc,mass,Q,xscri,coealphal,coealphah,icond2)
  if (isnan(epsi2)) then ! the equation becomes complex for values of rbar larger than rbaralpha: if calculated epsi2 is either complex or then negative, decrease the amunt by which icond2 is decreased below
    icond2 = icond2-0.01_iwp
    goto 644
  endif
!  print*, epsi1, epsi2
  649 continue
  nt = nt + 1
!  print*, nt, epsi3
  if (epsi1*epsi2.gt.0) then
  write(*,*) "Initial values for the bisection for inmstR not correct: "
  write(*,'(4(f10.4,a))') icond1, ' -> ', epsi1, ',  ', icond2, ' -> ', epsi2
  stop
  endif
  icond3 = (icond2+icond1)/2._iwp
  epsi3 = rbaralphaval(Kcmc,mass,Q,xscri,coealphal,coealphah,icond3)
  if (abs(epsi3).gt.1d-22) then !12
    if (epsi1*epsi3.gt.0) then
 icond1 = icond3
 epsi1 = epsi3
    else
 icond2 = icond3
 epsi2 = epsi3
    endif
    if (nt.gt.100000) stop "Arrived at 100000 iterations for rbaralpha - decrease accuracy."
    goto 649
  endif
  rbaralphacal = icond3
return
end function rbaralphacal
!!!! related function
function rbaralphaval(Kcmc,mass,Q,xscri,coealphal,coealphah,x)
implicit none
real(kind=iwp), intent(in)  :: Kcmc, mass, Q, xscri, coealphal, coealphah, x
real(kind=iwp) :: rbaralphaval, omega
omega = -Kcmc/6._iwp*(xscri+x)*(xscri-x)/xscri
rbaralphaval = coealphah + coealphal + Kcmc + (Kcmc**2*mass)/(9._iwp*Sqrt(-1 + (2*omega*mass)/x)) + (omega**2*mass)/(Sqrt(-1 + (2*omega*mass)/x)*x**2) - (2*omega*Sqrt(-1 + (2*omega*mass)/x))/x - (2*Kcmc**2*Sqrt(-1 + (2*omega*mass)/x)*x)/(9._iwp*omega) + (coealphah*Kcmc**2*x**2)/(9._iwp*omega**2) + (coealphal*Kcmc**2*x**2)/(9._iwp*omega**2) + (2*Kcmc**3*x**2)/(27._iwp*omega**2)
return
end function rbaralphaval

! not used, remove if not useful
  subroutine betacal(fac,x,dx,parameters,psic,aconfc)
  implicit none
  real(kind=iwp), dimension(:),   intent(in)  :: x
  type(parameters_type),intent(in)  :: parameters
  real(kind=wp), dimension(size(x)) :: psic
  real(kind=iwp) :: eaconf1, eaconf2, eaconf3, icond1, icond2, icond3
  real(kind=iwp) :: dx
  integer, parameter :: facnostag=2 ! do not change, to calculate inital data with nostag grid
  integer :: lb, li, ui, ub, nt, fac !, counter
  real(kind=iwp), dimension(:) :: aconfc !, allocatable

!    allocate(aconfc(1:size(x)*fac*facnostag))

    !Original number of points
    lb = lbound(x,1)
    ub = ubound(x,1)
    li = lb + parameters%grid%nghost
    ui = ub - parameters%grid%nghost
!print*, lb, li, ui, ub, dx
    !Increased number of points (to calculate with higher resolution)
    dx=dx/fac/facnostag
    ui = li + parameters%grid%ncells*fac*facnostag - 1
    if(parameters%grid%origin.eq.'nostag'.or.parameters%grid%origin.eq.'misman') ui = li + parameters%grid%ncells*fac*facnostag
    ub = ui + parameters%grid%nghost
!print*, lb, li, ui, ub, dx
!-!   if (parameters%grid%origin.eq.'nostag') stop "Calculation of aconf not prepared for a non-staggered grid."  !-!

  icond1 = -6._iwp
  nt = 0
  if(abs(parameters%grid%xmin).gt.1d-12) stop "Leftmost boundary not set at the imposed xmin = 0."
!  call betasli(icond1,lb,li,ui,ub,0._iwp,dx,eaconf1,aconfc,parameters)
!print*, icond1, eaconf1

    !Back to original number of points
    dx=dx*fac*facnostag
    ui = li + parameters%grid%ncells - 1
    if(parameters%grid%origin.eq.'nostag'.or.parameters%grid%origin.eq.'misman') ui = li + parameters%grid%ncells
    ub = ui + parameters%grid%nghost
!print*, lb, li, ui, ub, dx
    psic(li:ui) = aconfc(li+int((fac-1)/2._wp)+int(facnostag/2._wp):ui:fac*facnostag) ! stag and misman
	if(parameters%grid%origin.eq.'nostag') psic(li:ui) = aconfc(li:ui:fac*facnostag) !-! ! nostag
  return
  end subroutine betacal
  
!!! main subroutine when solving the cK slicing condition: determines rbaralpha, betabar and aconfalpha (compactification factor)
  subroutine cKslicing(fac,x,parameters,aconfc,acoomec,betarc)
  implicit none
  real(kind=iwp), dimension(:),   intent(in)  :: x
  real(kind=wp), dimension(:),   intent(inout)  :: aconfc, betarc, acoomec
  type(parameters_type), intent(in)  :: parameters
  real(kind=iwp) :: mass, Kcmc, Q, rbaralpha, coealphal, coealphah, xscri
  real(kind=iwp), dimension(:) , allocatable :: betabar, alphabar, betarbar, aconfbar, acoomebar, rcomp, aconf, betar, acoome!, alpha !, aconfbarlong, rcomplong
!  real(kind=wp), dimension(:) , allocatable :: aconfc
  integer :: fac, maxcount, counter, lb, ub, li, ui
  real(kind=iwp) :: res, err
  integer :: longer = 5
!  integer :: facrk4 = 2 ! it should depend on the runge-kutta4, so there's no need to make it easy to change the value
  
  maxcount = parameters%grid%ncells*fac+1 ! here calculating with a different grid, as things will have to be interpolated afterwards anyway

  allocate(betabar(1:2*maxcount)) !facrk4*
  allocate(alphabar(1:2*maxcount)) !facrk4*
  allocate(betarbar(1:2*maxcount)) !facrk4*
  
  mass = parameters%physics%mass
  Kcmc = parameters%physics%Kcmc
  Q = parameters%physics%Q
  coealphal = parameters%physics%coealphal
  coealphah = parameters%physics%coealphah
  xscri = parameters%physics%xscri

  rbaralpha = rbaralphacal(Kcmc,mass,Q,xscri,coealphal,coealphah)
  
  call betasli(2*maxcount,betabar,alphabar,betarbar,rbaralpha,parameters) !facrk4*
  
  allocate(rcomp(1:maxcount))
  allocate(aconfbar(1:maxcount))
  allocate(acoomebar(1:maxcount))
  
  call rcompcoord(maxcount,alphabar,rbaralpha,rcomp,aconfbar,acoomebar,parameters) !(1:facrk4*maxcount:facrk4)
  ! is the correct value of alpha set when calculating rcomp??? - check convergence!!! 

    lb = lbound(x,1)
    ub = ubound(x,1)
    li = lb + parameters%grid%nghost
    ui = ub - parameters%grid%nghost
    
if (.false.) then ! apparently not needed for good interpolation (here I would add a few points beyond both extrema)
!  allocate(rcomplong(1-longer:maxcount+longer))
!  allocate(aconfbarlong(1-longer:maxcount+longer))
 !   rcomplong(1:maxcount) = rcomp
  !  if (parameters%grid%origin.eq."nostag") rcomplong(1-longer:0) = -rcomp(longer+1:2:-1)
   ! if (parameters%grid%origin.eq."nostag") rcomplong(maxcount+1:maxcount+longer) = xscri + (xscri - rcomp(maxcount-1:maxcount-longer:-1))
    !aconfbarlong(1:maxcount) = aconfbar
!    if (parameters%grid%origin.eq."nostag") aconfbarlong(1-longer:0) = -aconfbar(longer+1:2:-1)
 !   if (parameters%grid%origin.eq."nostag") aconfbarlong(maxcount+1:maxcount+longer) = -aconfbar(maxcount-1:maxcount-longer:-1)
  !  do counter = 1-longer, maxcount+longer
!    print*, rcomplong(counter), aconfbarlong(counter)
   ! write(856,*) rcomplong(counter), aconfbarlong(counter)
    !enddo
endif

  allocate(aconf(lb:ub))
  allocate(betar(lb:ub))
  allocate(acoome(lb:ub))
!  allocate(alpha(lb:ub))

  do counter = li, ui
!  call NevilleCD(rcomplong, aconfbarlong, maxcount+2*longer, x(counter), res, err)
  call NevilleCD(rcomp, aconfbar, maxcount, x(counter), aconf(counter), err)
!  call NevilleCD(rcomp, alphabar, maxcount, x(counter), alpha(counter), err) ! alpha has to be calculated in a different way!
  call NevilleCD(rcomp, betarbar, maxcount, x(counter), betar(counter), err)
  call NevilleCD(rcomp, acoomebar, maxcount, x(counter), acoome(counter), err)
  aconfc(counter) = aconf(counter)  
  betarc(counter) = betar(counter)  
  acoomec(counter) = acoome(counter)  
!print*, x(counter), aconf(counter), aconfc(counter), beta(counter), betac(counter) !, alpha(counter)
write(955,*) x(counter), aconf(counter), aconfc(counter), betar(counter), betarc(counter), acoome(counter), acoomec(counter)  !, alpha(counter)
  enddo
  return
  end subroutine cKslicing  

!!!!!!!!!!!!!! calculation of beta (betabar) on the rbar equispaced grid, extending from rbaralpha to xscri !!!!! 
  subroutine betasli(maxcount,beta,alpha,betar,rbaralpha,parameters) 
  implicit none
  type(parameters_type), intent(in)  :: parameters
  real(kind=iwp) :: eaconf, icond
  real(kind=iwp) :: mass, Kcmc, Q, rbaralpha, coealphal, coealphah, xscri
  real(kind=iwp) :: dx, xv
  integer :: counter, maxcount
  real(kind=iwp) :: betap, betam, betarhs, k1, k2, k3, k4
  real(kind=iwp), dimension(1:maxcount) :: beta, alpha, betar

  mass = parameters%physics%mass
  Kcmc = parameters%physics%Kcmc
  Q = parameters%physics%Q
  coealphal = parameters%physics%coealphal
  coealphah = parameters%physics%coealphah
  xscri = parameters%physics%xscri

  dx = (xscri-rbaralpha)/real(maxcount-1,iwp)
  xv = rbaralpha
  icond = -(Kcmc*(-xv**2 + xscri**2)*Sqrt(-1 - (Kcmc*mass*(-xv**2 + xscri**2))/(3._iwp*xv*xscri)))/(6._iwp*xscri)
  betap = icond
  beta(1) = icond
  alpha(1) = sqrt(icond**2+(-Kcmc/6._iwp*(xscri+xv)*(xscri-xv)/xscri)**2*(1-2*mass*(-Kcmc/6._iwp*(xscri+xv)*(xscri-xv)/xscri)/xv))
  betar(1) = beta(1)*alpha(1)/( (-Kcmc/6._iwp*(xscri+xv)*(xscri-xv)/xscri) - xv*(Kcmc/3._iwp*xv/xscri) )
!print*, xv, betap, alpha(1)
write(755,*) xv, betap, alpha(1), betar(1)
  do counter = 2, maxcount
  betam = betap
  call betaprime(xv,betam,mass,Kcmc,Q,xscri,coealphal,coealphah,betarhs)
  if(.true.)then
!RK 4th order
  k1 = dx*betarhs
  xv = xv+dx/2._iwp
  call betaprime(xv,betam+k1/2._iwp,mass,Kcmc,Q,xscri,coealphal,coealphah,betarhs)
  k2 = dx*betarhs
  call betaprime(xv,betam+k2/2._iwp,mass,Kcmc,Q,xscri,coealphal,coealphah,betarhs)
  k3 = dx*betarhs
  xv = xv+dx/2._iwp
  call betaprime(xv,betam+k3,mass,Kcmc,Q,xscri,coealphal,coealphah,betarhs)
  k4 = dx*betarhs
  betap = betam + (k1+2*k2+2*k3+k4)/6._iwp
  else
!1st order
  betap = betam + dx*betarhs
  xv = xv + dx
  endif
  beta(counter) = betap
  alpha(counter) = sqrt(betap**2+(-Kcmc/6._iwp*(xscri+xv)*(xscri-xv)/xscri)**2*(1-2*mass*(-Kcmc/6._iwp*(xscri+xv)*(xscri-xv)/xscri)/xv))
  betar(counter) = beta(counter)*alpha(counter)/( (-Kcmc/6._iwp*(xscri+xv)*(xscri-xv)/xscri) - xv*(Kcmc/3._iwp*xv/xscri) )
!print*, xv, betap, alpha(counter)
write(755,*) xv, betap, alpha(counter), betar(counter) ! apparently ok
  if (abs(betap-betam).gt.1) betap = betam ! keep??
  end do
  return
  end subroutine betasli
!!!! related function
  subroutine betaprime(x,beta,mass,Kcmc,Q,xscri,coealphal,coealphah,betarhs)
  implicit none
  real(kind=iwp) :: x, mass, Kcmc, Q, xscri, coealphal, coealphah
  real(kind=iwp) :: betarhs, omega, domega, ddomega, beta
  omega = -Kcmc/6._iwp*(xscri+x)*(xscri-x)/xscri
  domega = Kcmc/3._iwp*x/xscri
  ddomega = Kcmc/3._iwp/xscri
  betarhs = (-coealphal - Kcmc - (2*omega*coealphah*mass)/x - (coealphah*Kcmc**2*x**2)/(9._iwp*omega**2) -    (coealphal*Kcmc**2*x**2)/(9._iwp*omega**2) - (2*Kcmc**3*x**2)/(27._iwp*omega**2) - (omega*mass*beta)/x**2 + (2*beta)/x +    (2*Kcmc**2*x*beta)/(9._iwp*omega**2) - (domega*beta)/(omega - domega*x) -    (domega*Kcmc**2*x**2*beta)/(9._iwp*omega**2*(omega - domega*x)) + (coealphah*beta**2)/omega**2 +    (domega*beta**3)/(omega**2*(omega - domega*x)) +    (coealphal*Sqrt(9 + (Kcmc**2*x**2)/omega**2)*Sqrt(1 - (2*omega*mass)/x + beta**2/omega**2))/3._iwp)/ (-(omega/(omega - domega*x)) - (Kcmc**2*x**2)/(9._iwp*omega*(omega - domega*x)) + beta**2/(omega*(omega - domega*x)))
  if (abs(xscri-x).lt.1d-7) betarhs = ((-9*beta**2*coealphah)/(-9*beta**2 + Kcmc**2*x**2) +  (18*beta**3*ddomega)/(domega*(-9*beta**2 + Kcmc**2*x**2)) + (9*beta**3)/(x*(-9*beta**2 + Kcmc**2*x**2)) -  (18*beta**2*coealphah*ddomega*x)/(domega*(-9*beta**2 + Kcmc**2*x**2)) - (3*beta*Kcmc**2*x)/(-9*beta**2 + Kcmc**2*x**2) +  (2*coealphah*Kcmc**2*x**2)/(-9*beta**2 + Kcmc**2*x**2) + (2*coealphal*Kcmc**2*x**2)/(-9*beta**2 + Kcmc**2*x**2) +  (4*Kcmc**3*x**2)/(3._iwp*(-9*beta**2 + Kcmc**2*x**2)) -  (3*beta**2*coealphal*Kcmc**2*x**2)/(Sqrt(beta**2*Kcmc**2*x**2)*(-9*beta**2 + Kcmc**2*x**2)) -  (beta*Kcmc**2*(4*domega**2*x + 6*ddomega*domega*x**2))/(domega**2*(-9*beta**2 + Kcmc**2*x**2)) +  (3*coealphah*Kcmc**2*x**2 + 3*coealphal*Kcmc**2*x**2 + 2*Kcmc**3*x**2 - 9*coealphal*Sqrt(beta**2*Kcmc**2*x**2))/  (3._iwp*(-9*beta**2 + Kcmc**2*x**2)) + (2*ddomega*x*    (3*coealphah*Kcmc**2*x**2 + 3*coealphal*Kcmc**2*x**2 + 2*Kcmc**3*x**2 - 9*coealphal*Sqrt(beta**2*Kcmc**2*x**2)))/  (3._iwp*domega*(-9*beta**2 + Kcmc**2*x**2)))/ (-1 - (27*beta**2)/(-9*beta**2 + Kcmc**2*x**2) + (18*beta*coealphah*x)/(-9*beta**2 + Kcmc**2*x**2) +  (3*Kcmc**2*x**2)/(-9*beta**2 + Kcmc**2*x**2) + (3*beta*coealphal*Kcmc**2*x**3)/(Sqrt(beta**2*Kcmc**2*x**2)*(-9*beta**2 + Kcmc**2*x**2))) ! rhs adapted to be evaluated at scri - omegas properly removed and limit calculated
  return
  end subroutine betaprime

!!!!!!!!!!!!!! calculation of rcomp (the compactified radial coordinate) and aconfalpha (the compactification factor) on the rbar equispaced grid !!!!! 
  subroutine rcompcoord(maxcount,alpha,rbaralpha,rcomp,aconf,acoome,parameters) 
  implicit none
  type(parameters_type), intent(in)  :: parameters
  real(kind=iwp) :: eaconf, icond
  real(kind=iwp) :: Kcmc, xscri, rbaralpha
  real(kind=iwp) :: dx, xv
  integer :: counter, maxcount
  real(kind=iwp) :: rcompp, rcompm, rrhs, k1, k2, k3, k4
  real(kind=iwp), dimension(1:maxcount) :: rcomp, aconf, acoome!, rbar
  real(kind=iwp), dimension(1:) :: alpha

  Kcmc = parameters%physics%Kcmc
  xscri = parameters%physics%xscri

  dx = -(xscri-rbaralpha)/real(maxcount-1,iwp)
  xv = xscri
  icond = 1._iwp
  rcompp = icond
!  rbar(maxcount) = xv
  rcomp(maxcount) = icond
  aconf(maxcount) = 0._iwp
  acoome(maxcount) = 1._iwp
!print*, xv, rcompp, aconf(maxcount)
write(855,*) xv, rcompp, aconf(maxcount), acoome(maxcount), alpha(maxcount)
  do counter = maxcount-1, 2, -1
  rcompm = rcompp
  call rprime(xv,rcompm,Kcmc,xscri,alpha(2*counter+1),rrhs) !!!index???
  if(.true.)then
!RK 4th order
  k1 = dx*rrhs
  xv = xv+dx/2._iwp
  call rprime(xv,rcompm+k1/2._iwp,Kcmc,xscri,alpha(2*counter+0),rrhs)!!!index???
  k2 = dx*rrhs
  call rprime(xv,rcompm+k2/2._iwp,Kcmc,xscri,alpha(2*counter+0),rrhs)!!!index???
  k3 = dx*rrhs
  xv = xv+dx/2._iwp
  call rprime(xv,rcompm+k3,Kcmc,xscri,alpha(2*counter-1),rrhs)!!!index???
  k4 = dx*rrhs
  rcompp = rcompm + (k1+2*k2+2*k3+k4)/6._iwp
  else
!1st order
  rcompp = rcompm + dx*rrhs
  xv = xv + dx
  endif
!  rbar(counter) = xv
  rcomp(counter) = rcompp
  aconf(counter) = rcomp(counter)/xv*(-Kcmc/6._iwp*(xscri+xv)*(xscri-xv)/xscri)
  acoome(counter) = rcomp(counter)/xv
!print*, xv, rcompp, aconf(counter)
write(855,*) xv, rcompp, aconf(counter), acoome(counter), alpha(2*counter-1)
  if (abs(rcompp-rcompm).gt.1) rcompp = rcompm ! keep??
  end do
!  rbar(1) = xv+dx
  rcomp(1) = 0._iwp
  aconf(1) = 0._iwp
  acoome(1) = 0._iwp
!print*, xv+dx, rcomp(1), aconf(1)
write(855,*) xv+dx, rcomp(1), aconf(1), acoome(1), alpha(1)
  return
  end subroutine rcompcoord
!!!! related function
  subroutine rprime(x,rcomp,Kcmc,xscri,alpha,rrhs)
  implicit none
  real(kind=iwp) :: x, Kcmc, xscri, alpha
  real(kind=iwp) :: rrhs, omega, domega, rcomp
  omega = -Kcmc/6._iwp*(xscri+x)*(xscri-x)/xscri
  domega = Kcmc/3._iwp*x/xscri
  rrhs = (omega-x*domega)/x/alpha*rcomp
  return
  end subroutine rprime

!! interpolation routine using the Neville algorithm, original copied from https://www.experts-exchange.com/questions/21197849/FORTRAN-Neville's-algorithm.html
   SUBROUTINE Neville(A, B, xm, xt, res, err)
     REAL(kind=iwp), dimension(:), intent( in ) :: A, B
     REAL(kind=iwp), intent( in ) :: xt
     REAL(kind=iwp), intent( out ) :: res, err
     INTEGER, intent ( in ) :: xm
     REAL(kind=iwp), dimension(xm) :: P
     INTEGER j,k
     P(:) = B(:)
       DO j=1, xm
          DO k=1, xm - j
               P(k)=((xt - A(k+j)) * P(k) + (A(k)-xt) * P(k+1))/(A(k)-A(k+j))
            END DO
       END DO
       res=P(1)
       err=P(1)-P(2)
       return
   END SUBROUTINE Neville
!! improved version following Numerical Recipes p. 119, then fixed using p. 103 from Numerical Recipes for FORTRAN
     SUBROUTINE NevilleCD(xa, ya, n, x, y, dy)
     INTEGER, intent ( in ) :: n
     REAL(kind=iwp), dimension(1:n), intent( in ) :: xa, ya
     REAL(kind=iwp), intent( in ) :: x
     REAL(kind=iwp), intent( out ) :: y, dy
     integer :: order 
     REAL(kind=iwp), dimension(1:n) :: c, d
     INTEGER m, i, ns
     real(kind=iwp) :: ho, hp, w, den
     
     order = 20 !n-1 !5 ! for a reason I don't understand I cannot set order=n-1 (interpolation using all points in the grid), because bad  values appear near the origin and at scri. 
     ns = minloc(abs(xa-x),1)
       c = ya
       d = ya
       
     y = ya(ns)
     ns = ns - 1
       DO m=1, order-1 
          DO i=1, n-m 
               ho = xa(i) - x
               hp = xa(i+m) - x
               w = c(i+1) - d(i)
               den = ho-hp
               if (abs(den).lt.1d-10) stop "Two points in rcomp grid are the same."
               den = w/den
               d(i) = hp*den
               c(i) = ho*den
            END DO
            if ((2*ns).lt.(n-m)) then 
              dy = c(ns+1)
            else
              dy = d(ns)
              ns = ns - 1
            endif
            y = y + dy
!            print*, m, x, res0, y
       END DO
       return
   END SUBROUTINE NevilleCD
!! improved version following Numerical Recipes p. 119, then fixed using p. 103 from Numerical Recipes for FORTRAN   
     SUBROUTINE NevilleCDdouble(xa, ya, n, x, y, dy)
     INTEGER, intent ( in ) :: n
     REAL(kind=wp), dimension(1:n), intent( in ) :: xa, ya
     REAL(kind=wp), intent( in ) :: x
     REAL(kind=wp), intent( out ) :: y, dy
     integer :: order 
     REAL(kind=wp), dimension(1:n) :: c, d
     INTEGER m, i, ns
     real(kind=wp) :: ho, hp, w, den, dif, dift
     
     order = 4 !6 !20 works fing for the interior, but not close to the boundaries !5 ! for a reason I don't understand I cannot set order=n-1 (interpolation using all points in the grid), because bad  values appear near the origin and at scri. 
     ns = minloc(abs(x-xa),1) 
!The previous line gives the same value of ns than the following do loop.
!     print*, "first", ns
!     ns = 1
!     dif = abs(x-xa(1))
!     do i = 1, n
!     	dift = abs(x-xa(i))
!     	if (dift.lt.dif) then
!     		ns = i
!     		dif = dift
!     	endif
!     enddo
!     print*, "second", ns
       c = ya
       d = ya
       
     y = ya(ns)
     ns = ns - 1
       DO m=1, order-1 
          DO i=1, n-m 
               ho = xa(i) - x
               hp = xa(i+m) - x
               w = c(i+1) - d(i)
               den = ho-hp
               if (abs(den).lt.1d-10) stop "Two points in rcomp grid are the same."
               den = w/den
               d(i) = hp*den
               c(i) = ho*den
            END DO
            if ((2*ns).lt.(n-m)) then 
              dy = c(ns+1)
            else
              dy = d(ns)
              ns = ns - 1
            endif
            y = y + dy
       END DO
       return
   END SUBROUTINE NevilleCDdouble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Cordoba - extreme Reissner-Nordstrom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ernrho(icond,x,mass)
  implicit none
  real(kind=wp) :: icond, x, mass, ernrho
  ernrho = x-(icond-mass**2/icond+2*mass*log(icond/mass))
  return
  end function ernrho

  function dernrho(icond,x,mass)
  implicit none
  real(kind=wp) :: icond, x, mass, dernrho
  dernrho = -(1+mass/icond)**2
  return
  end function dernrho
!!! Cordoba - extreme Reissner-Nordstrom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!! Initial data for Klein-Gordon field on the hyperboloidal foliation with a black hole !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine epsiomtbh(icond,lb,li,ui,ub,x,dx,erpsi,out,parameters,omegai,psicompl,psiA)
!  use numservicef90
  implicit none
  real(kind=iwp) :: erpsi, icond
  real(kind=iwp), intent(in)  :: x
  type(parameters_type), intent(in)  :: parameters
  real(kind=iwp) :: rhs, dx, a, center, sigma, xv, phikg, dphikg, pikg, mAr, ao, sigmao, centero, aim, centerim, sigmaim, imphikg, dimphikg, impikg
  real(kind=iwp) :: mass, Kcmc, Q, Ccmc, qp, masskgf, inmstR, aa, xscri, rmatchmax
  real(kind=iwp), dimension(2) :: psi, psim, prhs, k1, k2, k3, k4
  integer :: lb, ub, li, ui, counter, counterpsic, kgf, fac
  real(kind=iwp), dimension(lb:ub):: out
  real(kind=iwp), dimension(:), intent(in) :: psicompl, psiA, omegai
  mass = parameters%physics%mass
  Kcmc = parameters%physics%Kcmc
  Ccmc = parameters%physics%Ccmc
  Q = parameters%physics%Q
  qp = parameters%physics%qp
  masskgf = parameters%physics%masskgf
  aa = parameters%physics%aa
  xscri = parameters%physics%xscri
  rmatchmax = parameters%slice%rmatchmax
  if (aa.lt.0) aa = -3/Kcmc

  inmstR = inmstRcal(Kcmc,mass,Q,Ccmc)

  kgf = parameters%physics%kgf
  xv = x

  a = parameters%id%a; center = parameters%id%center; sigma = parameters%id%sigma
  aim = parameters%id%aim; centerim = parameters%id%centerim; sigmaim = parameters%id%sigmaim
  ao = parameters%id%ao; centero = parameters%id%centero; sigmao = parameters%id%sigmao

  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
!  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Kcmc,Ccmc,omegai((li-2)*2),psicompl((li-2)*2),parameters%physics%propsign)
  pikg = 0
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = 0
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
!  rhs = (24*dphikg**2*icond**9*inmstR**4*kgf*pi)/(63*Ccmc**2 + 9*Ccmc**2*icond**8 - 48*icond**8*inmstR**4 - icond**8*inmstR**6*Kcmc**2 + 5*icond**12*inmstR**6*Kcmc**2 )!- 8*icond**8*inmstR**2*Sqrt(9*Ccmc**2 + 9*inmstR**4 + 6*Ccmc*inmstR**3*Kcmc + inmstR**6*Kcmc**2 - 18*inmstR**3*mass) ! is supposed to exactly zero
!  rhs = (24*dphikg**2*icond**9*inmstR**4*kgf*pi)/(9*Ccmc**2*icond**8 - 48*icond**8*inmstR**4 - icond**8*inmstR**6*Kcmc**2 + 5*icond**12*inmstR**6*Kcmc**2  + 63*Ccmc**2*psiA(1+parameters%grid%nghost)**2)!for Arr*psiA  !-  8*icond**8*inmstR**2*Sqrt(9*Ccmc**2 + 9*inmstR**4 + 6*Ccmc*inmstR**3*Kcmc + inmstR**6*Kcmc**2 - 18*inmstR**3*mass)
!not generalized for aa and xscri  rhs = (-2*icond*inmstR**3*(-192*dphikg**2*icond**8*inmstR*kgf*Pi + 24*Ccmc*Kcmc*psiA(1+parameters%grid%nghost) + 2*inmstR**3*Kcmc**2*psiA(1+parameters%grid%nghost)**2))/(1008*Ccmc**2 + 144*Ccmc**2*icond**8 - 768*icond**8*inmstR**4 - 16*icond**8*inmstR**6*Kcmc**2 + 80*icond**12*inmstR**6*Kcmc**2 + 168*Ccmc*inmstR**3*Kcmc*psiA(1+parameters%grid%nghost) + 7*inmstR**6*Kcmc**2*psiA(1+parameters%grid%nghost)**2) ! for Arr+psiA ! - 128*icond**8*inmstR**2*Sqrt(9*Ccmc**2 + 9*inmstR**4 + 6*Ccmc*inmstR**3*Kcmc + inmstR**6*Kcmc**2 - 18*inmstR**3*mass)
  rhs = (-6*icond*inmstR**3*(-24*aa*Ccmc*psiA(1+parameters%grid%nghost) - 64*aa**2*dphikg**2*icond**8*inmstR*kgf*Pi*xscri + 6*inmstR**3*psiA(1+parameters%grid%nghost)**2*xscri))/(xscri*(1008*aa**2*Ccmc**2 + 144*aa**2*Ccmc**2*icond**8 - 768*aa**2*icond**8*inmstR**4 - 16*aa**2*icond**8*inmstR**6*Kcmc**2 + 80*aa**2*icond**12*inmstR**6*Kcmc**2 - 504*aa*Ccmc*inmstR**3*psiA(1+parameters%grid%nghost)*xscri + 63*inmstR**6*psiA(1+parameters%grid%nghost)**2*xscri**2))
  if ((abs(mass).le.1d-12).and.(abs(Ccmc).le.1d-12))  rhs = -(dphikg**2*icond*kgf*Pi)/3._iwp - (kgf*Pi*pikg**2)/(3._iwp*icond**7) - psiA(1+parameters%grid%nghost)**2/(16._iwp*icond**7) - icond/xscri**2 + icond**5/xscri**2
!print*, rhs
  psi(1) = icond+xv*xv/2._iwp*(rhs)
  psi(2) = xv*(rhs) ! not correct, rhs also has to be derived!!!
  out(li) = psi(1)
!print*, xv, phikg
write(558,*) xv, psi(1), li, out(li)
!  xv = xv + dx !no

  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai((li+1-2)*2),psicompl((li+1-2)*2),parameters%physics%propsign)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai((li+1-2)*2),psicompl((li+1-2)*2),parameters%physics%propsign)

  do counter = li+1, ui+1 ! when changed also take care of counter in psicompl!!
  psim = psi
!print*, xv, psim
  counterpsic = (counter-2)*2 ! seems ok
!  print*, counter, prhs, psim, psicompl(counterpsic), psiA(counterpsic)
  call p2rhsomtbh(xv,psim,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic),psicompl(counterpsic),psiA(counterpsic-2),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*,
!print*, prhs
!print*, counter, counterpsic, psicompl(counterpsic)
!if(counter.eq.(li+1)) print*, counter, counterpsic, psicompl(counterpsic)
  if(.true.)then
!RK 4th order
  k1 = dx*prhs
  xv = xv+dx/2._iwp
  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai(counterpsic+1),psicompl(counterpsic+1),parameters%physics%propsign)
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai(counterpsic+1),psicompl(counterpsic+1),parameters%physics%propsign)
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  call p2rhsomtbh(xv,psim+k1/2._iwp,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic+1),psicompl(counterpsic+1),psiA(counterpsic-1),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*, prhs
  k2 = dx*prhs
  call p2rhsomtbh(xv,psim+k2/2._iwp,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic+1),psicompl(counterpsic+1),psiA(counterpsic-1),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*, prhs
  k3 = dx*prhs
  xv = xv+dx/2._iwp
  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  call p2rhsomtbh(xv,psim+k3,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic+2),psicompl(counterpsic+2),psiA(counterpsic+0),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*, prhs
  k4 = dx*prhs
  psi = psim + kgf*(k1+2*k2+2*k3+k4)/6._iwp
  else
!1st order
  psi = psim + kgf*dx*prhs
  xv = xv + dx
  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  endif
!  print*, counter, prhs, psi!, psicompl(counterpsic+2)
  out(counter) = psi(1)
write(558,*) xv, psi(1), counter, out(counter), psicompl(counterpsic+2), prhs
!if (counter.eq.ui-1) print*, xv-dx, psi(1), prhs
!print*, counter, xv!, phikg
!  if (abs(psi(1)-psim(1)).gt.1) psi = psim ! new - good enough, NaNs do not appear anymore
  if (isnan(psi(1))) psi = psim ! new - good enough, NaNs do not appear anymore
  !if (abs(psi(1)).gt.1d100) psi = psim ! test 2026/06/22
!  if (abs(psi(1)-psim(1)).gt.1) then
!	psi = psim ! new - good enough, NaNs do not appear anymore
!	erpsi = psi(1) -1
!	return
!  endif
!print*, counter, psi(1), psim(1)
  end do
  erpsi = (psi(1)+psim(1))/2._iwp-1 !psi(1)**4 - 1+4*pi/3._iwp*dphikg*dphikg !psi(2) ! condition ! fix
!  write(*,*) psi, erpsi
!print*,
write(558,*)
  return
  end subroutine epsiomtbh

  subroutine epsiomtnostagbh(icond,lb,li,ui,ub,x,dx,erpsi,out,parameters,omegai,psicompl,psiA)
!  use numservicef90
  implicit none
  real(kind=iwp) :: erpsi, icond
  real(kind=iwp), intent(in)  :: x
  type(parameters_type), intent(in)  :: parameters
  real(kind=iwp) :: rhs, dx, a, center, sigma, xv, phikg, dphikg, pikg, mAr, ao, centero, sigmao, aim, centerim, sigmaim, imphikg, dimphikg, impikg
  real(kind=iwp) :: mass, Kcmc, Q, Ccmc, qp, masskgf, inmstR, aa, xscri, rmatchmax
  real(kind=iwp), dimension(2) :: psi, psim, prhs, k1, k2, k3, k4
  integer :: lb, ub, li, ui, counter, counterpsic, kgf, fac
  real(kind=iwp), dimension(lb:ub):: out
  real(kind=iwp), dimension(:), intent(in) :: psicompl, psiA, omegai
  mass = parameters%physics%mass
  Kcmc = parameters%physics%Kcmc
  Ccmc = parameters%physics%Ccmc
  Q = parameters%physics%Q
  qp = parameters%physics%qp
  masskgf = parameters%physics%masskgf
  aa = parameters%physics%aa
  xscri = parameters%physics%xscri
  rmatchmax = parameters%slice%rmatchmax
  if (aa.lt.0) aa = -3/Kcmc

  inmstR = inmstRcal(Kcmc,mass,Q,Ccmc)

  kgf = parameters%physics%kgf
  xv = x

  a = parameters%id%a; center = parameters%id%center; sigma = parameters%id%sigma
  aim = parameters%id%aim; centerim = parameters%id%centerim; sigmaim = parameters%id%sigmaim
  ao = parameters%id%ao; centero = parameters%id%centero; sigmao = parameters%id%sigmao

  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = 0
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = 0
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  rhs = (-6*icond*inmstR**3*(-24*aa*Ccmc*psiA(1+parameters%grid%nghost) - 64*aa**2*dphikg**2*icond**8*inmstR*kgf*Pi*xscri + 6*inmstR**3*psiA(1+parameters%grid%nghost)**2*xscri))/(xscri*(1008*aa**2*Ccmc**2 + 144*aa**2*Ccmc**2*icond**8 - 768*aa**2*icond**8*inmstR**4 - 16*aa**2*icond**8*inmstR**6*Kcmc**2 + 80*aa**2*icond**12*inmstR**6*Kcmc**2 - 504*aa*Ccmc*inmstR**3*psiA(1+parameters%grid%nghost)*xscri + 63*inmstR**6*psiA(1+parameters%grid%nghost)**2*xscri**2))
  if ((abs(mass).le.1d-12).and.(abs(Ccmc).le.1d-12))  rhs = -(dphikg**2*icond*kgf*Pi)/3._iwp - (kgf*Pi*pikg**2)/(3._iwp*icond**7) - psiA(1+parameters%grid%nghost)**2/(16._iwp*icond**7) - icond/xscri**2 + icond**5/xscri**2
!print*, rhs
  psi(1) = icond+xv*xv/2._iwp*(rhs)
  psi(2) = xv*(rhs) ! not correct, rhs also has to be derived!!!
  out(li+1) = psi(1) ! fixed for nostag
!print*, xv, phikg
write(558,*) xv, psi(1), li+1, out(li+1)
!  xv = xv + dx !no

  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai((li+1-2)*2),psicompl((li+1-2)*2),parameters%physics%propsign)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai((li+1-2)*2),psicompl((li+1-2)*2),parameters%physics%propsign)

  do counter = li+2, ui ! stag li+1, ui+1 ! when changed also take care of counter in psicompl!!
  psim = psi
!print*, xv, psim
  counterpsic = (counter-2)*2 -1 ! change seems ok
!  print*, counter, prhs, psim, psicompl(counterpsic), psiA(counterpsic)
  call p2rhsomtbh(xv,psim,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic),psicompl(counterpsic),psiA(counterpsic-2),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*,
!print*, prhs
!print*, counter, counterpsic, psicompl(counterpsic)
!if(counter.eq.(li+1)) print*, counter, counterpsic, psicompl(counterpsic)
  if(.true.)then
!RK 4th order
  k1 = dx*prhs
  xv = xv+dx/2._iwp
  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai(counterpsic+1),psicompl(counterpsic+1),parameters%physics%propsign)
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai(counterpsic+1),psicompl(counterpsic+1),parameters%physics%propsign)
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  call p2rhsomtbh(xv,psim+k1/2._iwp,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic+1),psicompl(counterpsic+1),psiA(counterpsic-1),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*, prhs
  k2 = dx*prhs
  call p2rhsomtbh(xv,psim+k2/2._iwp,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic+1),psicompl(counterpsic+1),psiA(counterpsic-1),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*, prhs
  k3 = dx*prhs
  xv = xv+dx/2._iwp
  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  call p2rhsomtbh(xv,psim+k3,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic+2),psicompl(counterpsic+2),psiA(counterpsic+0),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*, prhs
  k4 = dx*prhs
  psi = psim + kgf*(k1+2*k2+2*k3+k4)/6._iwp
  else
!1st order
  psi = psim + kgf*dx*prhs
  xv = xv + dx
  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  endif
!  print*, counter, prhs, psi!, psicompl(counterpsic+2)
  out(counter) = psi(1)
write(558,*) xv, psi(1), counter, out(counter), psicompl(counterpsic+2), prhs
  if (isnan(psi(1))) psi = psim
  end do
  erpsi = psi(1)-1
write(558,*)
  return
  end subroutine epsiomtnostagbh
  
    subroutine epsiomtmismanbh(icond,lb,li,ui,ub,x,dx,erpsi,out,parameters,omegai,psicompl,psiA)
!  use numservicef90
  implicit none
  real(kind=iwp) :: erpsi, icond
  real(kind=iwp), intent(in)  :: x
  type(parameters_type), intent(in)  :: parameters
  real(kind=iwp) :: rhs, dx, a, center, sigma, xv, phikg, dphikg, pikg, mAr, ao, centero, sigmao, aim, centerim, sigmaim, imphikg, dimphikg, impikg
  real(kind=iwp) :: mass, Kcmc, Q, Ccmc, qp, masskgf, inmstR, aa, xscri, rmatchmax
  real(kind=iwp), dimension(2) :: psi, psim, prhs, k1, k2, k3, k4
  integer :: lb, ub, li, ui, counter, counterpsic, kgf, fac
  real(kind=iwp), dimension(lb:ub):: out
  real(kind=iwp), dimension(:), intent(in) :: psicompl, psiA, omegai
  mass = parameters%physics%mass
  Kcmc = parameters%physics%Kcmc
  Ccmc = parameters%physics%Ccmc
  Q = parameters%physics%Q
  qp = parameters%physics%qp
  masskgf = parameters%physics%masskgf
  aa = parameters%physics%aa
  xscri = parameters%physics%xscri
  rmatchmax = parameters%slice%rmatchmax
  if (aa.lt.0) aa = -3/Kcmc

  inmstR = inmstRcal(Kcmc,mass,Q,Ccmc)

  kgf = parameters%physics%kgf
  xv = x

  a = parameters%id%a; center = parameters%id%center; sigma = parameters%id%sigma
  aim = parameters%id%aim; centerim = parameters%id%centerim; sigmaim = parameters%id%sigmaim
  ao = parameters%id%ao; centero = parameters%id%centero; sigmao = parameters%id%sigmao

  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = 0
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = 0
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  rhs = (-6*icond*inmstR**3*(-24*aa*Ccmc*psiA(1+parameters%grid%nghost) - 64*aa**2*dphikg**2*icond**8*inmstR*kgf*Pi*xscri + 6*inmstR**3*psiA(1+parameters%grid%nghost)**2*xscri))/(xscri*(1008*aa**2*Ccmc**2 + 144*aa**2*Ccmc**2*icond**8 - 768*aa**2*icond**8*inmstR**4 - 16*aa**2*icond**8*inmstR**6*Kcmc**2 + 80*aa**2*icond**12*inmstR**6*Kcmc**2 - 504*aa*Ccmc*inmstR**3*psiA(1+parameters%grid%nghost)*xscri + 63*inmstR**6*psiA(1+parameters%grid%nghost)**2*xscri**2))
  if ((abs(mass).le.1d-12).and.(abs(Ccmc).le.1d-12))  rhs = -(dphikg**2*icond*kgf*Pi)/3._iwp - (kgf*Pi*pikg**2)/(3._iwp*icond**7) - psiA(1+parameters%grid%nghost)**2/(16._iwp*icond**7) - icond/xscri**2 + icond**5/xscri**2 !! generalise for Q
!print*, rhs
  psi(1) = icond+xv*xv/2._iwp*(rhs)
  psi(2) = xv*(rhs) ! not correct, rhs also has to be derived!!!
  out(li) = psi(1) 
!print*, xv, phikg
write(558,*) xv, psi(1), li, out(li)
!  xv = xv + dx !no

  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai((li+1-2)*2),psicompl((li+1-2)*2),parameters%physics%propsign)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai((li+1-2)*2),psicompl((li+1-2)*2),parameters%physics%propsign)

  do counter = li+1, ui !??? ! nostag li+2, ui ! stag li+1, ui+1 ! when changed also take care of counter in psicompl!!
  psim = psi
!print*, xv, psim
  counterpsic = (counter-2)*2 !-1 ! change seems ok
!  print*, counter, prhs, psim, psicompl(counterpsic), psiA(counterpsic)
  call p2rhsomtbh(xv,psim,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic),psicompl(counterpsic),psiA(counterpsic-2),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*,
!print*, prhs
!print*, counter, counterpsic, psicompl(counterpsic)
!if(counter.eq.(li+1)) print*, counter, counterpsic, psicompl(counterpsic)
  if(.true.)then
!RK 4th order
  k1 = dx*prhs
  xv = xv+dx/2._iwp
  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai(counterpsic+1),psicompl(counterpsic+1),parameters%physics%propsign)
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai(counterpsic+1),psicompl(counterpsic+1),parameters%physics%propsign)
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  call p2rhsomtbh(xv,psim+k1/2._iwp,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic+1),psicompl(counterpsic+1),psiA(counterpsic-1),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*, prhs
  k2 = dx*prhs
  call p2rhsomtbh(xv,psim+k2/2._iwp,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic+1),psicompl(counterpsic+1),psiA(counterpsic-1),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*, prhs
  k3 = dx*prhs
  xv = xv+dx/2._iwp
  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  call p2rhsomtbh(xv,psim+k3,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omegai(counterpsic+2),psicompl(counterpsic+2),psiA(counterpsic+0),phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
!print*, prhs
  k4 = dx*prhs
  psi = psim + kgf*(k1+2*k2+2*k3+k4)/6._iwp
  else
!1st order
  psi = psim + kgf*dx*prhs
  xv = xv + dx
  phikg = kgf*valuephiomtbh(xv,a,center,sigma)
  dphikg = kgf*valuedphiomtbh(xv,a,center,sigma)
  pikg = kgf*valuepiomtbh(xv,a,center,sigma,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  imphikg = kgf*valueimphiomtbh(xv,aim,centerim,sigmaim)
  dimphikg = kgf*valuedimphiomtbh(xv,aim,centerim,sigmaim)
  impikg = kgf*valueimpiomtbh(xv,aim,centerim,sigmaim,mass,Q,Kcmc,Ccmc,omegai(counterpsic+2),psicompl(counterpsic+2),parameters%physics%propsign)
  mAr = parameters%physics%em*valuemAr(xv,ao,centero,sigmao)
  endif
!  print*, counter, prhs, psi!, psicompl(counterpsic+2)
  out(counter) = psi(1)
write(558,*) xv, psi(1), counter, out(counter), psicompl(counterpsic+2), prhs
  if (isnan(psi(1))) psi = psim
  end do
  erpsi = psi(1)-1
write(558,*)
  return
  end subroutine epsiomtmismanbh

  function valuephiomtbh(x,a,center,sigma)
  implicit none
  real(kind=iwp) :: x, valuephiomtbh, a, center, sigma
  valuephiomtbh =  (a*exp(-(x**2-center**2)**2/(4._iwp*sigma**4)))
  return
  end function valuephiomtbh

  function valuedphiomtbh(x,a,center,sigma)
  implicit none
  real(kind=iwp) :: x, valuedphiomtbh, a, center, sigma
  valuedphiomtbh = (-a*exp(-(x**2-center**2)**2/(4._iwp*sigma**4))/sigma**4*x*(x**2-center**2))
  return
  end function valuedphiomtbh

  function valuepiomtbh(x,a,center,sigma,mass,Q,Kcmc,Ccmc,omega,psic,propsign)
  implicit none
  real(kind=wp) :: propsign
  real(kind=iwp) :: x, valuepiomtbh, a, center, sigma, mass, Q, Kcmc, Ccmc, psic, omega
  valuepiomtbh = propsign*valuedphiomtbh(x,a,center,sigma) *( (3*Ccmc*psic**3+Kcmc*x**3)*omega**2 / ( psic**2 *sqrt( 9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*Ccmc*Kcmc*psic**3*x**3 - 18*mass*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6 ) ) ) !Schw: 9*Ccmc**2*psic**6+6*(Ccmc*Kcmc-3*mass)*psic**3*x**3+9*psic**2*x**4+Kcmc**2*x**6 
  return
  end function valuepiomtbh
  
  function valueimphiomtbh(x,a,center,sigma)
  implicit none
  real(kind=iwp) :: x, valueimphiomtbh, a, center, sigma
  valueimphiomtbh =  (a*exp(-(x**2-center**2)**2/(4._iwp*sigma**4)))
  return
  end function valueimphiomtbh

  function valuedimphiomtbh(x,a,center,sigma)
  implicit none
  real(kind=iwp) :: x, valuedimphiomtbh, a, center, sigma
  valuedimphiomtbh = (-a*exp(-(x**2-center**2)**2/(4._iwp*sigma**4))/sigma**4*x*(x**2-center**2))
  return
  end function valuedimphiomtbh

  function valueimpiomtbh(x,a,center,sigma,mass,Q,Kcmc,Ccmc,omega,psic,propsign)
  implicit none
  real(kind=wp) :: propsign
  real(kind=iwp) :: x, valueimpiomtbh, a, center, sigma, mass, Q, Kcmc, Ccmc, psic, omega
  valueimpiomtbh = propsign*valuedphiomtbh(x,a,center,sigma) *( (3*Ccmc*psic**3+Kcmc*x**3)*omega**2 / ( psic**2 *sqrt( 9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*Ccmc*Kcmc*psic**3*x**3 - 18*mass*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6 ) ) ) 
  return
  end function valueimpiomtbh
  
  function valuemAr(x,a,center,sigma)
  implicit none
  real(kind=iwp) :: x, valuemAr, a, center, sigma
  valuemAr = (a*x*exp(-(x**2-center**2)**2/(4._iwp*sigma**4)))
  return
  end function valuemAr

  subroutine p2rhsomtbh(x,psim,mass,Q,Kcmc,Ccmc,qp,masskgf,aa,omega,psic,psiA,phikg,dphikg,pikg,imphikg,dimphikg,impikg,mAr,prhs)
  implicit none
  real(kind=iwp) :: x, mass, Q, Kcmc, Ccmc, qp, masskgf, phikg, dphikg, pikg, imphikg, dimphikg, impikg, rhs, psic, psiA, aa, omega, dmEr, mPhi, mAr, Vkgval
  real(kind=iwp), dimension(2) :: psim, prhs
!  print*, "Warning, solver for psi not generalized to consider non-vanishing mEr or mPhi."
  dmEr = 0._iwp !the main part proportional to Q is included by hand !! divided by psi^6 -- generalize expression!! 
  mPhi = 0._iwp !! divided by psi^6 -- generalize!! 
  rhs = - pi*psim(1)*(dphikg*dphikg) - pi*psic**4*(pikg*pikg)/psim(1)**7/omega**4
  prhs(1) = psim(2)
  prhs(2) = rhs + (-3*Ccmc**2*psic**4*psiA**2)/(4._iwp*x**6*psim(1)**7) + (- Kcmc**2/(12._iwp*psic**2) + (3*Ccmc**2*psic**4)/(4._iwp*x**6))*psim(1) + (Kcmc**2*psim(1)**5)/(12._iwp*psic**2) + (-(1/x) - Sqrt(9*Ccmc**2*psic**6 + 6*Ccmc*Kcmc*psic**3*x**3 - 18*mass*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6)/(3._iwp*psic*x**3))*psim(2)
!print*, prhs
  ! at scri (when aconf(=psic)is zero) (using l'Hopital rule)
!  if (abs(psic).lt.1d-6) prhs(2) = rhs + (-3*Ccmc**2*psic**4)/(4._wp*x**6*psim(1)**7) + (3*Ccmc**2*psic**4)/(4._wp*x**6)*psim(1) - (1/x)*psim(2) + other terms
 ! if (abs(psic).lt.1d-6) prhs(2) = (4*(2*dphikg**2*pi*x*psim(1) + 2*psim(2) - 15*x*psim(1)**3*psim(2)**2))/(3._wp*x*(-1 + 5*psim(1)**4)) ! psic set exactly to zero
  if (abs(psic).lt.1d-6) prhs(2) = (4*(2*dphikg**2*pi*x*psim(1)))/(3._iwp*x*(-1 + 5*psim(1)**4)) ! psic and psim(2) set exactly to zero

! original before generalizing for Q nonzero
!  prhs(2) = rhs + (-3*omega**2*psiA**2)/(16._iwp*psic**2*psim(1)**7) - (3*Ccmc**2*psic**4)/(4._iwp*x**6*psim(1)**7) + (3*omega*Ccmc*psiA*psic)/(4._iwp*x**3*psim(1)**7) + (-Kcmc**2/(12._iwp*psic**2) + (3*Ccmc**2*psic**4)/(4._iwp*x**6))*psim(1) + (Kcmc**2*psim(1)**5)/(12._iwp*psic**2) + (-(1/x) - Sqrt(9*Ccmc**2*psic**6 + 6*Ccmc*Kcmc*psic**3*x**3 - 18*mass*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6)/(3._iwp*psic*x**3))*psim(2)
!  if (abs(psic).lt.1d-6) prhs(2) = (27*psiA**2 + 16*aa**2*dphikg**2*Kcmc**2*pi*psim(1)**8)/(6._iwp*aa**2*Kcmc**2*psim(1)**7*(-1 + 5*psim(1)**4)) ! psic, omega and psim(2) set exactly to zero
!  print*, x, masskgf

! original with masskg:  prhs(2) = (-(omega**2*(3*psiA**2 + 16*mPhi**2*phikg**2*pi*qp**2)) + psic**6*((-16*pi*pikg**2)/omega**4 - (12*Ccmc**2)/x**6) + (12*Ccmc*omega*psiA*psic**3)/x**3)/(16._iwp*psic**2*psim(1)**7) - (psic**3*Q + dmEr*omega**3*x**2)**2/(4._iwp*psic**4*x**4*psim(1)**3) +   ((-12*dphikg**2*pi - Kcmc**2/psic**2 + (9*Ccmc**2*psic**4)/x**6 + (3*psic**2*Q**2)/x**4)*psim(1))/12._iwp + ((Kcmc**2 - 12*masskgf**2*phikg**2*pi)*psim(1)**5)/(12._iwp*psic**2) + (-(1/x) - Sqrt(9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*(Ccmc*Kcmc - 3*mass)*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6)/(3._iwp*psic*x**3))*psim(2) -((omega**4*phikg**2*pi*qp**2*mAr**2*psim(1)**9)/psic**4)
Vkgval = Vkgfval(x,masskgf,omega,phikg,imphikg)
!write(562,*) x, phikg, Vkgval
!  prhs(2) = (-3*omega**2*psiA**2)/(16._wp*psic**2*psim(1)**7) - (pi*pikg**2*psic**4)/(omega**4*psim(1)**7) - (mPhi**2*omega**2*phikg**2*pi*qp**2)/(psic**2*psim(1)**7) - (3*Ccmc**2*psic**4)/(4._wp*x**6*psim(1)**7) + (3*Ccmc*omega*psiA*psic)/(4._wp*x**3*psim(1)**7) -  (dmEr**2*omega**6)/(4._wp*psic**4*psim(1)**3) - (psic**2*Q**2)/(4._wp*x**4*psim(1)**3) - (dmEr*omega**3*Q)/(2._wp*psic*x**2*psim(1)**3) - dphikg**2*pi*psim(1) - (Kcmc**2*psim(1))/(12._wp*psic**2) + (3*Ccmc**2*psic**4*psim(1))/(4._wp*x**6) + (psic**2*Q**2*psim(1))/(4._wp*x**4) + (Kcmc**2*psim(1)**5)/(12._wp*psic**2) - (pi*Vkgval*psim(1)**5)/psic**2 - (omega**4*phikg**2*pi*qp**2*mAr**2*psim(1)**9)/psic**4 - psim(2)/x - (Sqrt(9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*(Ccmc*Kcmc - 3*mass)*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6)*psim(2))/(3._wp*psic*x**3) ! used until 2026/08/16
  prhs(2) = (-3*omega**2*psiA**2)/(16._iwp*psic**2*psim(1)**7) - (impikg**2*pi*psic**4)/(omega**4*psim(1)**7) - (pi*pikg**2*psic**4)/(omega**4*psim(1)**7) + (2*impikg*mPhi*phikg*pi*psic*qp)/(omega*psim(1)**7) - (2*imphikg*mPhi*pi*pikg*psic*qp)/(omega*psim(1)**7) - (imphikg**2*mPhi**2*omega**2*pi*qp**2)/(psic**2*psim(1)**7) - (mPhi**2*omega**2*phikg**2*pi*qp**2)/(psic**2*psim(1)**7) - (3*Ccmc**2*psic**4)/(4._iwp*x**6*psim(1)**7) + (3*Ccmc*omega*psiA*psic)/(4._iwp*x**3*psim(1)**7) - (dmEr**2*omega**6)/(4._iwp*psic**4*psim(1)**3) - (psic**2*Q**2)/(4._iwp*x**4*psim(1)**3) - (dmEr*omega**3*Q)/(2._iwp*psic*x**2*psim(1)**3) - dimphikg**2*pi*psim(1) - dphikg**2*pi*psim(1) - (Kcmc**2*psim(1))/(12._iwp*psic**2) + (3*Ccmc**2*psic**4*psim(1))/(4._iwp*x**6) + (psic**2*Q**2*psim(1))/(4._iwp*x**4) + (Kcmc**2*psim(1)**5)/(12._iwp*psic**2) + (2*dphikg*imphikg*omega**2*pi*qp*mAr*psim(1)**5)/psic**2 - (2*dimphikg*omega**2*phikg*pi*qp*mAr*psim(1)**5)/psic**2 - (imphikg**2*omega**4*pi*qp**2*mAr**2*psim(1)**9)/psic**4 - (omega**4*phikg**2*pi*qp**2*mAr**2*psim(1)**9)/psic**4 - psim(2)/x - (Sqrt(9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*(Ccmc*Kcmc - 3*mass)*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6)*psim(2))/(3._iwp*psic*x**3) - (pi*psim(1)**5*Vkgval)/psic**2
!  print*,1, impikg, dimphikg
!  print*,2, pikg, dphikg
  return
  end subroutine p2rhsomtbh

  subroutine p2rhsomtpsiA(x,psiA,mass,Q,Kcmc,Ccmc,qp,aa,omega,psic,phikg,dphikg,pikg,imphikg,dimphikg,impikg,prhs)
  implicit none
  real(kind=iwp) :: x, mass, Q, Kcmc, Ccmc, aa, phikg, dphikg, pikg, imphikg, dimphikg, impikg, rhs, psic, omega, qp, mAr, mPhi
  real(kind=iwp) :: psiA, prhs
  ! mPhi is assumed to be zero, as otherwise there would be a term depending on mAr
  mAr = 0
  mPhi = 0
  rhs = -(8*dphikg*pi*pikg*psic**3)/omega**3
!  prhs = rhs + psiA*(-(1/x) + Sqrt(omega**2*aa**4 + aa**2*x**2)/(omega*aa**2*x) - Sqrt(9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*Ccmc*Kcmc*psic**3*x**3 - 18*mass*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6)/(psic*x**3)) ! used until 2026/08/16
  !(-8*dphikg*pi*pikg*psic**3)/omega**3 - (Kcmc*psiA*x)/(3._iwp*omega) - (psiA*Sqrt(9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*(Ccmc*Kcmc - 3*mass)*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6))/(psic*x**3)
  !rhs + psiA*(-(1/x) + Sqrt(omega**2*aa**4 + aa**2*x**2)/(omega*aa**2*x) - Sqrt(9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*Ccmc*Kcmc*psic**3*x**3 - 18*mass*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6)/(psic*x**3)) ! original with Q added
  prhs = (-8*dimphikg*impikg*pi*psic**3)/omega**3 - (8*dphikg*pi*pikg*psic**3)/omega**3 - (Kcmc*psiA*x)/(3._iwp*omega) - (psiA*Sqrt(9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*(Ccmc*Kcmc - 3*mass)*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6))/(psic*x**3) !- (8*impikg*mAr*phikg*pi*psic*qp*psim(1)**4)/omega + (8*imphikg*mAr*pi*pikg*psic*qp*psim(1)**4)/omega - 8*dphikg*imphikg*mPhi*pi*qp*psim(1)**6 + 8*dimphikg*mPhi*phikg*pi*qp*psim(1)**6 + (8*imphikg**2*mAr*mPhi*omega**2*pi*qp**2*psim(1)**10)/psic**2 + (8*mAr*mPhi*omega**2*phikg**2*pi*qp**2*psim(1)**10)/psic**2
  if (abs(psic).lt.1d-6) prhs = 0
!print*, omega, psiA, rhs, prhs
  return
  end subroutine p2rhsomtpsiA

!!!!!!!!!!!!!!!!!!!!!!!!!! Initial data for Klein-Gordon field on the hyperboloidal foliation with a black hole !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!! Initial data for Klein-Gordon field on relaxed BH trumpet geometry !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function valuephiomtbhdouble(x,a,center,sigma)
  implicit none
  real(kind=wp) :: x, valuephiomtbhdouble, a, center, sigma
  valuephiomtbhdouble =  (a*exp(-(x**2-center**2)**2/(4._iwp*sigma**4)))
  return
  end function valuephiomtbhdouble

  function valuedphiomtbhdouble(x,a,center,sigma)
  implicit none
  real(kind=wp) :: x, valuedphiomtbhdouble, a, center, sigma
  valuedphiomtbhdouble = (-a*exp(-(x**2-center**2)**2/(4._iwp*sigma**4))/sigma**4*x*(x**2-center**2))
  return
  end function valuedphiomtbhdouble
  
  subroutine p2rhsrelaxedtrumpet(x,psim,Kcmc,xscri,z4,u,du,ddu,phikg,dphikg,prhs)
  implicit none
  real(kind=wp) :: x, Kcmc, xscri, phikg, dphikg, rhs, omega, domega, ddomega
  real(kind=wp), dimension(2) :: psim, prhs
  real(kind=wp), dimension(:) :: u, du, ddu
  real(kind=wp) :: chi, grr, Arr, K, Lambdar, alpha, betar, Theta, dchi, dgrr, dArr, dK, dLambdar, dalpha, dbetar, dTheta, ddchi, ddgrr, ddalpha, ddbetar
  integer :: z4
  omega = -Kcmc*(xscri+x)*(xscri-x)/(6._wp*xscri)
  domega = Kcmc*x/(3._wp*xscri)
  ddomega = Kcmc/(3._wp*xscri) 
  
  chi    	= u(1)
  grr   	= u(2)
  Arr   	= u(4)
  K   		= u(5)
  Lambdar   	= u(6)
  alpha   	= u(7)
  betar   	= u(8)
  Theta	= u(12)
  dchi    	= du(1)
  dgrr   	= du(2)
  dArr   	= du(4)
  dK   	= du(5)
  dLambdar   	= du(6)
  dalpha   	= du(7)
  dbetar   	= du(8)
  dTheta	= du(12)
  ddchi    	= ddu(1)
  ddgrr   	= ddu(2)
  ddalpha   	= ddu(7)
  ddbetar   	= ddu(8)
  
  rhs = (-3*Arr**2)/(16._wp*chi*grr*psim(1)**7) - (5*dchi**2*psim(1))/(16._wp*chi**2) + (ddchi*psim(1))/(4._wp*chi) + (ddomega*psim(1))/(2._wp*omega) - (dchi*domega*psim(1))/(4._wp*omega*chi) - (3*domega**2*psim(1))/(4._wp*omega**2) - (15*dgrr**2*psim(1))/(64._wp*grr**2) + (ddgrr*psim(1))/(8._wp*grr) - (dchi*dgrr*psim(1))/(4._wp*chi*grr) - (dgrr*domega*psim(1))/(2._wp*omega*grr) - omega**2*dphikg**2*pi*psim(1) - 2*omega*domega*dphikg*phikg*pi*psim(1) - domega**2*phikg**2*pi*psim(1) - psim(1)/(4._wp*x**2) + (grr**1.5_wp*psim(1))/(4._wp*x**2) + (dchi*psim(1))/(2._wp*chi*x) + (domega*psim(1))/(omega*x) + (5*dgrr*psim(1))/(8._wp*grr*x) + (grr*K**2*psim(1)**5)/(12._wp*omega**2*chi) + (grr*K*Kcmc*psim(1)**5)/(6._wp*omega**2*chi) + (grr*Kcmc**2*psim(1)**5)/(12._wp*omega**2*chi) + (grr*K*Theta*z4*psim(1)**5)/(3._wp*omega**2*chi) + (grr*Kcmc*Theta*z4*psim(1)**5)/(3._wp*omega**2*chi) + (grr*Theta**2*z4*psim(1)**5)/(3._wp*omega**2*chi) + (dchi*psim(2))/(2._wp*chi) + (domega*psim(2))/omega + (dgrr*psim(2))/grr - (2*psim(2))/x
!  print*, rhs
  prhs(1) = psim(2)
  prhs(2) = rhs 
!  print*, x, psim, prhs
  return
  end subroutine p2rhsrelaxedtrumpet

  subroutine epsirelaxtrumpet(icond,a,center,sigma,lb,li,ui,ub,xdoa,xdob,x,dx,erpsi,out,parameters,xvect,u)
  use numservicef90
  implicit none
  real(kind=wp) :: icond
  real(kind=wp) :: erpsi
  real(kind=wp) :: x
  type(parameters_type), intent(in)  :: parameters
  real(kind=wp) :: rhs, dx, a, center, sigma, phikg, dphikg !, pikg
  real(kind=wp) :: Kcmc, xscri, xdoa, xdob, omega, domega, ddomega, err
  real(kind=wp), dimension(2) :: psi, psim, prhs, k1, k2, k3, k4
  integer :: lb, ub, li, ui, counter, counterpsic, kgf, fac, z4, last=12, ivar
  real(kind=wp), dimension(lb:ub) :: out
  real(kind=wp), dimension(:,:), intent(in) :: u
  real(kind=wp), dimension(size(u,1),size(u,2)) :: du, ddu
  real(kind=wp), dimension((2*(size(u(1,:))-2*parameters%grid%nghost)+2*parameters%grid%nghost)+1) :: chi, grr, Arr, K, Lambdar, alpha, betar, Theta, dchi, dgrr, dArr, dK, dLambdar, dalpha, dbetar, dTheta, ddchi, ddgrr, ddalpha, ddbetar
  real(kind=wp), dimension(lb:ub) :: xvect
  real(kind=wp), dimension((2*(size(u(1,:))-2*parameters%grid%nghost)+2*parameters%grid%nghost)+1) :: xinterp, chiinterp, ddchiinterp
  real(kind=wp), dimension(size(u,1),(2*(size(u(1,:))-2*parameters%grid%nghost)+2*parameters%grid%nghost)+1) :: uinterp, duinterp1, dduinterp1, duinterp2, dduinterp2, duinterp, dduinterp
  Kcmc = parameters%physics%Kcmc
  xscri = parameters%physics%xscri
  z4 = parameters%physics%z4
  
  omega = -Kcmc*(xscri+x)*(xscri-x)/(6._wp*xscri)
  domega = Kcmc*x/(3._wp*xscri)
  ddomega = Kcmc/(3._wp*xscri) 
  
  du = u
  ddu = u
!  print*, xdoa, xdob
  if (parameters%moldef%deriv_method == 'c8') then
     call mdiff_c8(du(1:last,:), xdoa, xdob, parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
     call mdiff_c6(du(1:last,:), xdoa, xdob, parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call mdiff_c4(du(1:last,:), xdoa, xdob, parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
     call mdiff(du(1:last,:), xdoa, xdob, parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
  
  if (parameters%moldef%deriv_method == 'c8') then
     call m2diff_c8(ddu(1:last,:), xdoa, xdob, parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
     call m2diff_c6(ddu(1:last,:), xdoa, xdob, parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call m2diff_c4(ddu(1:last,:), xdoa, xdob, parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
     call m2diff(ddu(1:last,:), xdoa, xdob, parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
  
  do counter = lb, ub
  	write(785,*) xvect(counter), ddu(1,counter), ddu(2,counter), du(4,counter), du(5,counter), du(12,counter)
  enddo
stop
!print*, "Checking matrix bounds. Current run: ", lb, li, ui, ub, " Data bounds: ", lbound(u,2), ubound(u,2), lbound(ddu,2), ubound(ddu,2) ! looks alright, they all start at 1, li is the same for this and the previous run


if (parameters%grid%origin.ne.'stag') stop "Interpolation of data only implemented for staggered grid."

!  do counter = lb, ub
 ! 	write(788,*) xvect(counter), chi(counter), dchi(counter), ddchi(counter), counter
  !enddo

 do counter = 1, 2*parameters%grid%ncells + 2*parameters%grid%nghost + 1 ! may have to update for different types of grids
  	xinterp(counter) = xdoa+(counter+0)*dx/2._wp 
  	do ivar = 1, size(u,1)
  		call NevilleCDdouble(xvect, u(ivar,:), size(u,2), xinterp(counter), uinterp(ivar,counter), err)
  		call NevilleCDdouble(xvect, du(ivar,:), size(u,2), xinterp(counter), duinterp1(ivar,counter), err)
  		call NevilleCDdouble(xvect, ddu(ivar,:), size(u,2), xinterp(counter), dduinterp1(ivar,counter), err)
  	enddo
  	write(790,*) xinterp(counter), (uinterp(ivar,counter),ivar=1,1)!size(u,1))
  	write(791,*) xinterp(counter), (duinterp1(ivar,counter),ivar=1,1)!size(u,1))
  	write(792,*) xinterp(counter), (dduinterp1(ivar,counter),ivar=1,1)!size(u,1))
  enddo

  
  duinterp2 = uinterp
  dduinterp2 = uinterp
  
  
! print*, lbound(xinterp,1), ubound(xinterp,1), xinterp(1), xinterp(lbound(xinterp,1)), xinterp(ubound(xinterp,1)) 
    if (parameters%moldef%deriv_method == 'c8') then
     call mdiff_c8(duinterp2(1:last,:), xinterp(lbound(xinterp,1)), xinterp(ubound(xinterp,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
     call mdiff_c6(duinterp2(1:last,:), xinterp(lbound(xinterp,1)), xinterp(ubound(xinterp,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call mdiff_c4(duinterp2(1:last,:), xinterp(lbound(xinterp,1)), xinterp(ubound(xinterp,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
     call mdiff(duinterp2(1:last,:), xinterp(lbound(xinterp,1)), xinterp(ubound(xinterp,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
  
  if (parameters%moldef%deriv_method == 'c8') then
     call m2diff_c8(dduinterp2(1:last,:), xinterp(lbound(xinterp,1)), xinterp(ubound(xinterp,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
     call m2diff_c6(dduinterp2(1:last,:), xinterp(lbound(xinterp,1)), xinterp(ubound(xinterp,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call m2diff_c4(dduinterp2(1:last,:), xinterp(lbound(xinterp,1)), xinterp(ubound(xinterp,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
     call m2diff(dduinterp2(1:last,:), xinterp(lbound(xinterp,1)), xinterp(ubound(xinterp,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
  
  do counter = 1, 2*parameters%grid%ncells + 2*parameters%grid%nghost + 1 ! may have to update for different types of grids
  	write(793,*) xinterp(counter), (duinterp2(ivar,counter),ivar=1,1)!size(u,1))
  	write(794,*) xinterp(counter), (dduinterp2(ivar,counter),ivar=1,1)!size(u,1))
  enddo
  
  duinterp = duinterp1 ! these look better near the boundaries
  dduinterp = dduinterp1

  chi    	= uinterp(1,:)
  grr   	= uinterp(2,:)
  Arr   	= uinterp(4,:)
  K   		= uinterp(5,:)
  Lambdar   	= uinterp(6,:)
  alpha   	= uinterp(7,:)
  betar   	= uinterp(8,:)
  Theta	= uinterp(12,:)
  dchi    	= duinterp(1,:)
  dgrr   	= duinterp(2,:)
  dArr   	= duinterp(4,:)
  dK   	= duinterp(5,:)
  dLambdar   	= duinterp(6,:)
  dalpha   	= duinterp(7,:)
  dbetar   	= duinterp(8,:)
  dTheta	= duinterp(12,:)
  ddchi    	= dduinterp(1,:)
  ddgrr   	= dduinterp(2,:)
  ddalpha   	= dduinterp(7,:)
  ddbetar   	= dduinterp(8,:)

  kgf = parameters%physics%kgf

  phikg = kgf*valuephiomtbhdouble(x,a,center,sigma)
  dphikg = kgf*valuedphiomtbhdouble(x,a,center,sigma)
!  pikg = kgf*valuepisymmetric(x,a,center,sigma,betar)

  rhs = 0 !-(icond/xscri**2) - (dphikg**2*icond*Kcmc**2*Pi*xscri**2)/36. + (3*icond**5)/(xscri**2*chi(li)) - (3*Arr(li)**2)/(16.*icond**7*chi(li)) - (3*icond*dchi(li)**2)/(16.*chi(li)**2) + (icond*ddchi(li))/(4.*chi(li)) + (icond*ddgrr(li))/8. + (icond*dchi(li)*dgrr(li))/(4.*chi(li)) + (17*icond*dgrr(li)**2)/64. + (6*icond**5*K(li))/(Kcmc*xscri**2*chi(li)) + (3*icond**5*K(li)**2)/(Kcmc**2*xscri**2*chi(li)) + (12*icond**5*z4*Theta(li))/(Kcmc*xscri**2*chi(li)) + (12*icond**5*z4*K(li)*Theta(li))/(Kcmc**2*xscri**2*chi(li)) + (12*icond**5*z4*Theta(li)**2)/(Kcmc**2*xscri**2*chi(li)) ! update: special l'hopitalized rhs at the origin
  
  psi(1) = icond + x*x/2._wp*(rhs)
  psi(2) = x*(rhs) ! not correct, rhs also has to be derived!!!
  out(li) = psi(1)
  write(558,*) x, psi(1), li, out(li)

do counter = li+1, ui+1 ! when changed also take care of counter in psicompl!!
  psim = psi
  counterpsic = (counter-2)*2 ! seems ok ! check if this also works for u 
  call p2rhsrelaxedtrumpet(x,psim,Kcmc,xscri,z4,uinterp(:,counterpsic),duinterp(:,counterpsic),dduinterp(:,counterpsic),phikg,dphikg,prhs)
  !if(counter.eq.(li+1)) print*, counter, counterpsic, u(1,counterpsic), u(7,counterpsic)
if(.true.)then
!RK 4th order
  k1 = dx*prhs
  x = x+dx/2._wp
  phikg = kgf*valuephiomtbhdouble(x,a,center,sigma)
  dphikg = kgf*valuedphiomtbhdouble(x,a,center,sigma)
  call p2rhsrelaxedtrumpet(x,psim+k1/2._wp,Kcmc,xscri,z4,uinterp(:,counterpsic+1),duinterp(:,counterpsic+1),dduinterp(:,counterpsic+1),phikg,dphikg,prhs)
  k2 = dx*prhs
  call p2rhsrelaxedtrumpet(x,psim+k2/2._wp,Kcmc,xscri,z4,uinterp(:,counterpsic+1),duinterp(:,counterpsic+1),dduinterp(:,counterpsic+1),phikg,dphikg,prhs)
  k3 = dx*prhs
  x = x+dx/2._wp
  phikg = kgf*valuephiomtbhdouble(x,a,center,sigma)
  dphikg = kgf*valuedphiomtbhdouble(x,a,center,sigma)
  call p2rhsrelaxedtrumpet(x,psim+k3,Kcmc,xscri,z4,uinterp(:,counterpsic+2),duinterp(:,counterpsic+2),dduinterp(:,counterpsic+2),phikg,dphikg,prhs)
  k4 = dx*prhs
  psi = psim + kgf*(k1+2*k2+2*k3+k4)/6._wp
else
!1st order
  psi = psim + kgf*dx*prhs
  x = x + dx
  phikg = kgf*valuephiomtbhdouble(x,a,center,sigma)
  dphikg = kgf*valuedphiomtbhdouble(x,a,center,sigma)
endif
  out(counter) = psi(1)
  write(558,*) x, psi(1), counter, out(counter), prhs
  !if (counter.eq.ui-1) print*, x-dx, psi(1), prhs
  !  if (abs(psi(1)-psim(1)).gt.1) psi = psim ! new - good enough, NaNs do not appear anymore
  if (isnan(psi(1))) psi = psim ! new - good enough, NaNs do not appear anymore
end do
  erpsi = (psi(1)+psim(1))/2._wp-1 !psi(1)**4 - 1+4*pi/3._iwp*dphikg*dphikg !psi(2) ! condition ! fix
  write(558,*)
  return
  end subroutine epsirelaxtrumpet

!!!!!!!!!!!!!!!!!!!!!!!!!! Initial data for Klein-Gordon field on relaxed BH trumpet geometry !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mol_subroutines
