!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The heart of the module is the subroutine molstep
! it does the time evolution step using method of lines
!
!   statevect = statevect + dt*rhs
!
! The rhs is calculated by the subroutine calcrhs
! calcrhs is where the physical system comes in.
! It will be modified when equations change.
!
! Because of the abstract nature of mol,
! this module doesn't need modification for different systems.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mol
  use prec
  use parameters_mod
  use molparameters
  use mol_subroutines

  implicit none

    real(kind=wp), dimension(:,:,:), private,   allocatable, save  :: tmp
    real(kind=wp), parameter :: third=1.0_wp/3.0_wp, sixth=1.0_wp/6.0_wp
    integer  :: tmpsteps

contains

  subroutine mol_param_testing(parameters)
    implicit none
    type(parameters_type), intent(in) :: parameters

    if(parameters%timer%dt <= 0.0_wp) then
  stop 'dt smaller than zero!'
    endif

    if(parameters%timer%t > parameters%timer%end_time) then
  stop 'end_time too small!'
    endif
end subroutine mol_param_testing

  subroutine molstep(statevect, constr, gridvect, gauge, parameters, psic)

    !use boundaries
    implicit none

    interface
  subroutine calcrhs(statevect, constr, gridvect, gauge, rhs, parameters, psic)
use parameters_mod

real(kind=wp), dimension(:,:),intent(in) :: statevect, gauge, constr
real(kind=wp), dimension(:),intent(in) :: gridvect, psic
type(parameters_type), intent(in) :: parameters

real(kind=wp), dimension(size(statevect,1), &
    &size(statevect,2)), intent(out) :: rhs
  end subroutine calcrhs
    end interface


    real(kind=wp), dimension(:,:),intent(inout) :: statevect, gauge, constr
    real(kind=wp), dimension(:),intent(inout) :: gridvect, psic
!real(kind=wp), dimension(:) ::  psic

    type(parameters_type), intent(inout) :: parameters

    real(kind=wp), dimension(size(statevect,1), size(statevect,2)) :: rhs

    integer  :: counter = 1 !used in the Adams-Bashford-Moulton method
    integer  :: modul  !used in the Adams-Bashford-Moulton method


    ! executable statements


    select case(parameters%moldef%int_method)

    case ('rk3')
  tmpsteps = 2
    case ('rk4')
  tmpsteps = 2
    case ('rk4_3')
  tmpsteps = 3
    end select


    if (.NOT.(allocated(tmp)) .AND.  (tmpsteps > 0))  then
  allocate(tmp(size(statevect,1), size(statevect,2), tmpsteps))
  print*, 'allocated ', tmpsteps, ' temporary time levels for MoL'
    endif

    select case(parameters%moldef%int_method)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Heun's third-order formula Trefethen 1994 p75
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case ('rk3')
! step 1
  call calcrhs(statevect, constr, gridvect, gauge, rhs, parameters, psic)
  tmp(:,:,1) = statevect + third * parameters%timer%dt * rhs
  call setboundaries(tmp(:,:,1), gridvect, parameters)

  parameters%timer%t = parameters%timer%t + third*parameters%timer%dt

! step 2
  call calcrhs(tmp(:,:,1), constr, gridvect, gauge, rhs, parameters, psic)
  tmp(:,:,2) = statevect + 2._wp * third * parameters%timer%dt * rhs
  call setboundaries(tmp(:,:,2), gridvect, parameters)

  parameters%timer%t = parameters%timer%t + third*parameters%timer%dt

! step 3
  call calcrhs(tmp(:,:,2), constr, gridvect, gauge, rhs, parameters, psic)

  statevect = 0.25_wp*(statevect + 3._wp*(tmp(:,:,1)+ parameters%timer%dt * rhs))
  call setboundaries(statevect, gridvect, parameters)

  parameters%timer%t = parameters%timer%t + third*parameters%timer%dt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This seems to be a correct implementation of rk4 Trefethen 1994 p 76
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case ('rk4_3')  ! rk4 with _3_ timelevels

! step 1
  call calcrhs(statevect, constr, gridvect, gauge, rhs, parameters, psic)
  tmp(:,:, 1) = statevect + rhs * parameters%timer%dt * 0.5_wp
  call setboundaries(tmp(:,:,1), gridvect, parameters)
  parameters%timer%t = parameters%timer%t + parameters%timer%dt *0.5_wp

  tmp(:,:, 3) = rhs

! step 2
  call calcrhs(tmp(:,:, 1), constr, gridvect, gauge, rhs, parameters, psic)
  tmp(:,:, 2) = statevect + rhs * parameters%timer%dt * 0.5_wp
  call setboundaries( tmp(:,:, 2), gridvect, parameters)

  tmp(:,:, 3) =  tmp(:,:, 3) + 2.0_wp * rhs

! step 3
  call calcrhs(tmp(:,:, 2), constr, gridvect, gauge, rhs, parameters, psic)

  tmp(:,:, 1) = statevect  + rhs * parameters%timer%dt
  call setboundaries(tmp(:,:, 1), gridvect, parameters)
  parameters%timer%t = parameters%timer%t + parameters%timer%dt *0.5_wp

  tmp(:,:, 3) =  tmp(:,:, 3) + 2.0_wp * rhs

! step 4
  call calcrhs(tmp(:,:, 1), constr, gridvect, gauge, rhs, parameters, psic)

  tmp(:,:, 3) =  tmp(:,:, 3) + rhs

  statevect = statevect + parameters%timer%dt*sixth *  tmp(:,:, 3)
  call setboundaries(statevect, gridvect, parameters)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    RK4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case ('rk4')  ! rk4 with 2 timelevels

  stop 'rk4 is not correctly implemented.'

! step 1
  call calcrhs(statevect, constr, gridvect, gauge, rhs, parameters, psic)
  tmp(:,:, 1) = statevect + rhs * parameters%timer%dt * 0.5_wp
  call setboundaries(tmp(:,:,1), gridvect, parameters)
  parameters%timer%t = parameters%timer%t + parameters%timer%dt *0.5_wp

  statevect = statevect + (parameters%timer%dt *sixth) * rhs

! step 2

  call calcrhs(tmp(:,:, 1), constr, gridvect, gauge, rhs, parameters, psic)

  tmp(:,:, 2) = statevect  + rhs * parameters%timer%dt * 0.5_wp
  call setboundaries(statevect, gridvect, parameters)
  statevect = statevect + (parameters%timer%dt / 3.0_wp) * rhs

! step 3
  call calcrhs(tmp(:,:, 2), constr, gridvect, gauge, rhs, parameters, psic)

  tmp(:,:, 1) = statevect  + rhs * parameters%timer%dt
  call setboundaries(statevect, gridvect, parameters)
  parameters%timer%t = parameters%timer%t + parameters%timer%dt *0.5_wp

  statevect = statevect + (parameters%timer%dt / 3.0_wp) * rhs

! step 4
  call calcrhs(tmp(:,:, 1), constr, gridvect, gauge, rhs, parameters, psic)

  statevect = statevect  + parameters%timer%dt / (6.0_wp) * rhs
  call setboundaries(statevect, gridvect, parameters)


    case DEFAULT
  stop 'method not supported'

    end select
  end subroutine molstep



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setting initial data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialdata(u, constr, xv, gauge, parameters, psic)
!  use parameters_mod
  use numservicef90
  implicit none

  real(kind=wp), dimension(:,:), intent(inout) :: u, gauge
  real(kind=wp), dimension(:,:),intent(out)    :: constr
  real(kind=iwp), dimension(:),   intent(inout)  :: xv
  real(kind=wp), dimension(size(xv)) :: rh
  real(kind=wp), dimension(size(xv)) :: x
  type(parameters_type),intent(inout)  :: parameters

  character :: cc*5,ac*5

  real(kind=wp) :: a, sigma, center, ao, centero, sigmao, aim, sigmaim, centerim !, v, h, k
  integer :: i, j, n, noise!, kgf


  real(kind=wp), dimension(size(u,1), size(u,2)) :: rhs
!
  real(kind=wp), dimension(size(u,1), size(u,2)) :: du, ddu, u_diss, dru, druin
  real(kind=wp), dimension(size(x)) :: chi, grr, gtt, Arr, K, Lambdar, alpha, betar, Br, derchi, dergrr, dergtt, deralpha, derbetar, phikg, pikg, Theta, Zr, Lambdazr, RR, TT, rephi, repi, imphi, impi, mEr, mAr, mPhi, mPsi
  real(kind=wp), dimension(size(x)) :: dchi, dgrr, dgtt, dArr, dK, dLambdar, dalpha, dbetar, dBr, dderchi, ddergrr, ddergtt, dderalpha, dderbetar, dphikg, dpikg, dTheta, dZr, dLambdazr, dRR, dTT, drephi, drepi, dimphi, dimpi, dmEr, dmAr, dmPhi, dmPsi
  real(kind=wp), dimension(size(x)) :: drchi, drgrr, drgtt, drArr, drK, drLambdar, dralpha, drbetar, drBr, drphikg, drpikg, drTheta, drZr, drRR, drTT, drrephi, drrepi, drimphi, drimpi, drmEr, drmAr, drmPhi, drmPsi ! for advection terms
  real(kind=wp), dimension(size(x)) :: ddchi, ddgrr, ddgtt, ddArr, ddK, ddLambdar, ddalpha, ddbetar, ddBr, ddphikg, ddpikg, ddTheta, ddZr, ddrephi, ddrepi, ddimphi, ddimpi, ddmEr, ddmAr, ddmPhi, ddmPsi
  real(kind=wp), dimension(size(x)) :: auxmass, auxp, auxpp, div, squ, omega, domega, ddomega, dddomega, root, root2, coefom, ders, der, fwe, fwea, lcoef, Vkgfunc, Vkgprimefunc  ! sqar, sqar2, sqar3, gsqarba, gsqarba2, ugsqarb, usqarb, gsqarb,
  real(kind=wp), dimension(size(x)) :: psic, dpsic, ddpsic, coefpsi, rootp, rootp2, beta
  real(kind=wp), dimension(size(x)) :: dersdalpha, dersdgtt, dersdchi, dersArr,dersK, dersLambdar, dersalphau, dedersgttchi, ddersgttchi, alphahat, dalphahat, betarhat, dbetarhat
  real(kind=wp), dimension(size(x)) :: grrb, dgrrb, ddgrrb, gttb, dgttb, ddgttb, coefalpha, dcoefalpha, ddcoefalpha, coefgrr, dcoefgrr, ddcoefgrr, coefgrrb, dcoefgrrb, ddcoefgrrb, nuplog, nconst
  real(kind=wp) :: v, eta, lambda, mu, eps, mass, cosm, cosmk, aa, cuplog, charm, xmax, k1, k2, xscri, coealphah, coealphal, coebetar, coeLambdar, coeZr, pphi, trafK, xexc, nuplogv, nconstv, masskgf, qp
  real(kind=wp) :: Kcmc, Ccmc, Q, propsign, rmatchmin, rmatchmax, tvaldiag
  integer	:: evollapse, evolshift, last, nro, nr0, kgf, null, counter, numr, z4, z4c, omt, imatchmin, imatchmax, em, gr
  real(kind=iwp) :: dxv, funcikg, xiv
  real(kind=wp) :: dx
  real(kind=wp), dimension(size(x)) :: transition, intslicing, extslicing, intshift, extshift
!!
  real(kind=wp), dimension(size(x)) :: psi0, psi, dpsi, ddpsi, psirhs, psidiff, ppsi, psiA
  real(kind=wp) :: psim, ppsim, icond1, icond2, icond3, prhs, epsi1, epsi2, epsi3
  integer :: nt
  integer, parameter :: fac = 1 !9 !27 ! the higher the convergence order (6 or 8), the higher the multiple of 3 to be set !!!!now set to 1 for psic - improve??
  real(kind=wp), dimension(:), allocatable :: out ! made allocatable to avoid segmentation faults
  real(kind=iwp), dimension(:), allocatable :: psicompl, psiAcompl, outi, psii2, omegai
  real(kind=iwp), dimension(size(x)) :: psii
  real(kind=wp) :: dt = 0.000001, rms = 1, rms0 = 0
  integer :: lb, ub, li, ui

  integer, parameter :: n2v = 18 

  real(kind=iwp) :: ai, centeri, sigmai, Kcmci, Ccmci, Qi, massi, aai, xscrii, aoi, centeroi, sigmaoi, aimi, centerimi, sigmaimi, qpi
  real(kind=iwp) :: icond1i, icond2i, icond3i, epsi1i, epsi2i, epsi3i
  real(kind=iwp) :: k1i, k2i, k3i, k4i, dphikgi, pikgi, prhsi, imphikgi, dimphikgi, impikgi

  real(kind=wp) :: rbaralpha
  real(kind=iwp) :: coealphahi, coealphali
  real(kind=wp), dimension(size(u,1), 2*(size(u(1,:))-2*parameters%grid%nghost)+2*parameters%grid%nghost) :: udo !, dudo, ddudo ! for reading in data for the perturbed relaxed trumpet
  real(kind=wp), dimension(2*(size(u(1,:))-2*parameters%grid%nghost)+2*parameters%grid%nghost) :: xdo ! for reading in the gridvector for the perturbed relaxed trumpet

  allocate(out(1:size(x)*fac))
  allocate(outi(1:size(x)*fac))
  allocate(psicompl(1:size(x)*fac*2))
  allocate(psiAcompl(1:size(x)*fac*2))
  allocate(omegai(1:size(x)*fac*2))
  allocate(psii2(1:size(x)*fac*2))

  x = xv

    ! executable statements

    lb = lbound(u, 2)
    ub = ubound(u, 2)

    li = lb + parameters%grid%nghost
    ui = ub - parameters%grid%nghost
!!
!print*, lb, li, ui, ub
  cuplog = parameters%slice%cuplog
  charm = parameters%slice%charm
  evollapse = parameters%slice%evollapse
  evolshift = parameters%slice%evolshift
  v = parameters%physics%v
  eta = parameters%physics%eta
  lambda = parameters%physics%lambda
  mu = parameters%physics%mu
  aai = parameters%physics%aa
  massi = parameters%physics%mass
  cosm = parameters%physics%cosm
  gr = parameters%physics%gr
  kgf = parameters%physics%kgf
  em = parameters%physics%em
  omt = parameters%physics%omt
  cosmk = cosm/3._wp
  xmax = parameters%grid%xmax
  xscrii = parameters%physics%xscri
  xexc = parameters%physics%xexc
  z4 = parameters%physics%z4
  z4c = parameters%physics%z4c
  k1 = parameters%physics%k1
  k2 = parameters%physics%k2
  coealphahi = parameters%physics%coealphah
  coealphali = parameters%physics%coealphal
  nuplogv = parameters%physics%nuplogv
  nconstv = parameters%physics%nconstv
  coebetar = parameters%physics%coebetar
  coeLambdar = parameters%physics%coeLambdar
  coeZr = parameters%physics%coeZr
  pphi = parameters%physics%pphi
  trafK = parameters%physics%trafK
  Kcmci = parameters%physics%Kcmc
  Ccmci = parameters%physics%Ccmc
  Qi = parameters%physics%Q
  qpi = parameters%physics%qp
  propsign = parameters%physics%propsign
  rmatchmin = parameters%slice%rmatchmin
  rmatchmax = parameters%slice%rmatchmax
  tvaldiag = parameters%slice%tvaldiag

  ! value of Ccmc calculated from Kcmc and mass that makes the foliations be trumpets
  if (Ccmci.lt.0) Ccmci = ccmccal(Kcmci,massi,Qi)

  if (aai.lt.0) aai = -3/Kcmci

  mass = massi; Q = Qi; Kcmc = Kcmci; Ccmc = Ccmci; aa = aai; xscri = xscrii; qp = qpi
  
  coealphah = coealphahi; coealphal = coealphali

  parameters%physics%Ccmc = Ccmci !! dangerous?

!xscrii = inmstRcal(Kcmci,massi,Qi,Ccmci) ! correctness check
!print*, Kcmci, massi, Qi, Ccmci !, xscrii

  eps = parameters%moldef%dissipation_eps
!
dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / dble(parameters%grid%ncells)
dxv   = abs(parameters%grid%xmax - parameters%grid%xmin) / real(parameters%grid%ncells,iwp)
!if (parameters%grid%origin.eq.'stag') dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / (dble(parameters%grid%ncells+1)) !change for staggered grid  !RUINS CONVERGENCE!!!!!!!
if (parameters%grid%origin.eq.'misman') then
	dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / (dble(parameters%grid%ncells)+0.5_wp)
	dxv   = abs(parameters%grid%xmax - parameters%grid%xmin) / (real(parameters%grid%ncells,iwp)+0.5_iwp)
endif

  if (.true.) then
  omega = (xscri+x)*(xscri-x)/(2._wp*aa*xscri) !-Kcmc*(xscri+x)*(xscri-x)/(6._wp*xscri)
  domega = -x/(aa*xscri) !Kcmc*x/(3._wp*xscri)
  ddomega = -1/(aa*xscri) !Kcmc/(3._wp*xscri)!  dddomega = 0
  coefom = omega-x*domega
  root2 = 9*omega**6*Ccmc**2 + 6*omega**3*(Ccmc*Kcmc - 3*mass)*x**3 + 9*omega**2*x**4 + Kcmc**2*x**6
  root = sqrt(root2)
  else
  omega = (xscri+x)*(xscri-x)/(2._wp*aa*xscri)
  domega = -x/(aa*xscri)
  ddomega = -1/(aa*xscri)
  dddomega = 0
  coefom = omega-x*domega
  root2 = x*x+aa*aa*omega*omega
  root = sqrt(root2)
  endif

  do counter = lbound(psicompl,1), ubound(psicompl,1)
   	omegai(counter) = (xscrii+dxv/2._wp*(counter-parameters%grid%nghost-1))*(xscrii-dxv/2._iwp*(counter-parameters%grid%nghost-1))/(2._iwp*aai*xscrii)
  enddo

if((abs(mass).gt.1d-12).and.(abs(Ccmc).gt.1d-12))then
  write(*,*) "Determination of the overall compactification factor:"
  call aconfcal(fac,xv,dxv,parameters,psic,psicompl)
	do counter = 1, abs(parameters%grid%nghost) ! bcs for aconf(=psic)
		psic(lbound(u, 2) + counter - 1) = - psic(lbound(u, 2) + 2*abs(parameters%grid%nghost) - counter)
		if (parameters%grid%origin.eq.'nostag') psic(lbound(u, 2) + counter - 1) = - psic(lbound(u, 2) + 2*abs(parameters%grid%nghost) - counter + 1)
		psic(ubound(u, 2)-counter+1) = omega(ubound(u, 2)-counter+1)
	enddo

  psicompl(ubound(psicompl, 1)-(parameters%grid%nghost*2-1+parameters%grid%nghost)+1) = (2*xscri+dx/2._wp)*(-dx/2._wp)/(2._wp*aa*xscri)
  psicompl(ubound(psicompl, 1)-(parameters%grid%nghost*2-1+parameters%grid%nghost)+2) = (2*xscri+dx)*(-dx)/(2._wp*aa*xscri)
  psicompl(ubound(psicompl, 1)-(parameters%grid%nghost*2-1+parameters%grid%nghost)+3) = (2*xscri+3*dx/2._wp)*(-3*dx/2._wp)/(2._wp*aa*xscri)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do counter = lbound(u, 2), ubound(u, 2)
	write(556,*) x(counter), psic(counter), omega(counter) !, psic(counter)/x(counter)*1.68423694_wp
  enddo
  do counter = lbound(psicompl, 1), ubound(psicompl, 1)
	write(557,*) dx/2._wp*(counter-parameters%grid%nghost-1), psicompl(counter)
  enddo
else
write(*,*) "As mass = 0 and Ccmc = 0, aconf will be set to omega."
  psic = omega
!  do counter = lbound(psicompl,1), ubound(psicompl,1)
 !  	psicompl(counter) = (xscrii+dx/2._wp*(counter-parameters%grid%nghost-1))*(xscrii-dx/2._wp*(counter-parameters%grid%nghost-1))/(2._wp*aai*xscrii)
  !enddo
  psicompl = omegai
  do counter = lbound(psicompl, 1), ubound(psicompl, 1)
	write(557,*) dx/2._wp*(counter-parameters%grid%nghost-1), psicompl(counter)
  enddo
endif

!  rootp2 = 9*Ccmc**2*psic**6+6*(Ccmc*Kcmc-3*mass)*psic**3*x**3+9*psic**2*x**4+Kcmc**2*x**6!9*Ccmc**2 + 6*(Ccmc*Kcmc-3*mass)*x**3*psic**3 + 9*x**4*psic**4 + Kcmc**2*x**6*psic**6
do i=lbound(rootp2,1), ubound(rootp2,1) ! for some unknown reason the previous expression gives rootp2=0
  rootp2(i) = 9*Ccmc**2*psic(i)**6+6*(Ccmc*Kcmc-3*mass)*psic(i)**3*x(i)**3+9*psic(i)**2*x(i)**4+Kcmc**2*x(i)**6
enddo
  rootp = sqrt(abs(rootp2))

!print*, "rootp", psic(100), rootp2(100), sqrt(abs(9*Ccmc**2*psic(100)**6+6*(Ccmc*Kcmc-3*mass)*psic(100)**3*x(100)**3+9*psic(100)**2*x(100)**4+Kcmc**2*x(100)**6)), rootp(100)

if(.false.)then
  if (.true.) then
  grrb = 9*x**4*coefpsi**2/rootp2 !9/(9 + Kcmc**2*psic**4*x**2) !(9*coefpsi**2*psic**6*x**4)/rootp2 ! 9*coefom*coefom*x**4/root2
  gttb = 1
  else
  grrb = coefom*coefom*aa*aa/root2
  gttb = 1
  endif
  if (parameters%moldef%deriv_method == 'c8') then
call vdiff_c8(grrb, dgrrb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
call vdiff_c6(grrb, dgrrb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
call vdiff_c4(grrb, dgrrb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
call vdiff(grrb, dgrrb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
  if (parameters%moldef%deriv_method == 'c8') then
call v2diff_c8(grrb, ddgrrb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
call v2diff_c6(grrb, ddgrrb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
call v2diff_c4(grrb, ddgrrb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
call v2diff(grrb, ddgrrb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
  if (parameters%moldef%deriv_method == 'c8') then
call vdiff_c8(gttb, dgttb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
call vdiff_c6(gttb, dgttb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
call vdiff_c4(gttb, dgttb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
call vdiff(gttb, dgttb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
  if (parameters%moldef%deriv_method == 'c8') then
call v2diff_c8(gttb, ddgttb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
call v2diff_c6(gttb, ddgttb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
call v2diff_c4(gttb, ddgttb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
call v2diff(gttb, ddgttb, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
else
  grrb = 1 !coefom*coefom*aa*aa/root2 !grrb
  dgrrb = 0 !aa*aa*2*coefom/root2*(-x*ddomega-coefom*(x+aa*aa*omega*domega)/root2) !dgrrb
  ddgrrb = 0 !aa*aa*2/root2*((x*ddomega+2*coefom*(x+aa*aa*omega*domega)/root2)*(x*ddomega+coefom*(x+aa*aa*omega*domega)/root2)-coefom/root2*((ddomega+x*dddomega)*root2+coefom*(1+aa*aa*domega*domega+aa*aa*omega*ddomega)-x*ddomega*(x+aa*aa*omega*domega)-2*coefom*(x+aa*aa*omega*domega)**2/root2)) !ddgrrb
!! to avoid extra noise coming from the numerical derivatives or expressions
  gttb = 1
  dgttb = 0
  ddgttb = 0
endif


  a = parameters%id%a
  sigma = parameters%id%sigma
  center = parameters%id%center
  aim = parameters%id%aim
  sigmaim = parameters%id%sigmaim
  centerim = parameters%id%centerim
  ao = parameters%id%ao
  sigmao = parameters%id%sigmao
  centero = parameters%id%centero

  noise = parameters%id%noise

!  h = parameters%slice%h
!  k = parameters%slice%k
!  v = parameters%physics%v
!  kgf = parameters%physics%kgf

  ! 1...pi
  ! 2...psi
  ! 3...phi

!alpha: + dble(1-kgf)*(a*exp(-(x**2-center**2)**2/(4._wp*sigma**4))) !evenperturb!
!dalpha: + dble(1-kgf)*(-a*exp(-(x**2-center**2)**2/(4._wp*sigma**4))/sigma**4*x*(x**2-center**2)) !evenperturb!
!phi: kgf*(a*exp(-(x**2-center**2)**2/(4._wp*sigma**4))) !evenperturb!

!do i = lb, ub
!write(*,*) x(i)
!enddo

  select case(parameters%id%id_type)

  case ('from_data')

open(13, FILE=trim(file_prefix)//'last_data.txt', STATUS='OLD')
read(13,*) n, parameters%timer%t
parameters%timer%end_time = parameters%timer%end_time + parameters%timer%t
if(n.ne.size(u(1,:))) then !parameters%grid%ncells+1
write(*,*) n
stop "wrong number of data points"
end if
write(*,*) parameters%grid%ncells+1-n,size(u(1,:))-n
!do i= 1, parameters%grid%ncells+1
do i = lbound(u,2), ubound(u,2)
read(13,'(19(e35.24e3))') x(i), (u(j,i),j=lbound(u,1), ubound(u,1)) 
enddo
close(13)

!!!!!!! statiobetar  
  include '../source.f90' ! to calculate the constraints 
!!!!!!! statiobetar

   case ('bh-hypcomp')

!linear!  	if (parameters%moldef%eqsystem.eq.'nonlinear') then

call random_number(u)
u = (u - 0.5_wp)*1.0d-8 !0.1_wp
!u(:,lbound(u,2)+parameters%grid%nghost:lbound(u,2)+parameters%grid%nghost+200) = 0.0_wp
!u(:,ubound(u,2)-parameters%grid%nghost-200:ubound(u,2)-parameters%grid%nghost) = 0.0_wp

!Kcmc=0._wp; omega=1._wp; domega=0._wp; ddomega=0._wp; psic=1._wp; dpsic=0._wp; ddpsic=0._wp


u(1,:) = noise*u(1,:) + psic**2/omega**2 
u(2,:) = noise*u(2,:) + grrb  
u(3,:) = noise*u(3,:) + (-2*Ccmc*psic**3)/(omega*x**3)  
u(4,:) = noise*u(4,:) + 0   
u(5,:) = noise*u(5,:) + 0 
u(6,:) = noise*u(6,:) + sqrt(abs( 9*Ccmc**2*psic**6 + 9*psic**4*Q**2*x**2 + 6*Ccmc*Kcmc*psic**3*x**3 - 18*mass*psic**3*x**3 + 9*psic**2*x**4 + Kcmc**2*x**6 ))*omega/(3._wp*psic*x**2) + dble(1-kgf)*(a*exp(-(x**2-center**2)**2/(4._wp*sigma**4)))
!original for Schw: sqrt(abs(9*Ccmc**2*psic**6+6*(Ccmc*Kcmc-3*mass)*psic**3*x**3+9*psic**2*x**4+Kcmc**2*x**6))*omega/(3._wp*psic*x**2) + dble(1-kgf)*(a*exp(-(x**2-center**2)**2/(4._wp*sigma**4)))
u(7,:) = noise*u(7,:) + ( Ccmc*psic**3/x**2 + Kcmc*x/3._wp )
u(8,:) = noise*u(8,:) + 0.0_wp

!Set to zero here, calculation of initial data if kgf is 1 done below
u(9,:) = noise*u(9,:) + kgf*(a*exp(-(x**2-center**2)**2/(4._wp*sigma**4)))/omega ! duplicated from below -- if gr evolved, it gets overwritten there ! real parts of scalar field
u(10,:) = noise*u(10,:)  -(-a*exp(-(x**2-center**2)**2/(4._wp*sigma**4))*x*(x**2-center**2)/sigma**4) * ( Ccmc*psic**3/x**2 + Kcmc*x/3._iwp )*(-propsign-1)/omega ! duplicated from below -- if gr evolved, it gets overwritten there ! real parts of scalar field
u(11,:) = noise*u(11,:) + kgf*(aim*exp(-(x**2-centerim**2)**2/(4._wp*sigmaim**4)))/omega ! duplicated from below -- if gr evolved, it gets overwritten there ! imaginary parts of scalar field
u(12,:) = noise*u(12,:) -(-aim*exp(-(x**2-centerim**2)**2/(4._wp*sigmaim**4))*x*(x**2-centerim**2)/sigmaim**4) * ( Ccmc*psic**3/x**2 + Kcmc*x/3._iwp )*(-propsign-1)/omega ! duplicated from below -- if gr evolved, it gets overwritten there ! imaginary parts of scalar field

u(13,:) = 0._wp + em*(ao*x*exp(-(x**2-centero**2)**2/(4._wp*sigmao**4))) ! Ar
u(14,:) = Q*u(1,:)**(1.5_wp)/x**2 ! Er
u(15,:) = 0._wp ! Phi
u(16,:) = 0._wp ! Psi

! diag variables for creating the conformal diagrams
u(17,:) = ( atan(tvaldiag + (x/omega) + Sqrt(9/Kcmc**2 + x**2/omega**2) + 3/Kcmc) - atan(tvaldiag -(x/omega) + Sqrt(9/Kcmc**2 + x**2/omega**2) + 3/Kcmc) )!/2._wp ! Minkowski
u(18,:) = ( atan(tvaldiag + (x/omega) + Sqrt(9/Kcmc**2 + x**2/omega**2) + 3/Kcmc) + atan(tvaldiag -(x/omega) + Sqrt(9/Kcmc**2 + x**2/omega**2) + 3/Kcmc) )!/2._wp  !Minkowski

if (gr.eq.1) then ! solving for the perturbation (not needed for Cowling)

!! Initial data KGF hypcomp with BH
! Setting initial data for chi using psi with a shooting&matching method
if (kgf.eq.1) then
write(*,*) "Initial data:"
write(*,*) "Warning: intial solver considers perturbations in both real and imaginary parts of the scalar field, but its electric charge is assumed to be zero in part of the solving procedure." !"Warning: only the real part of the scalar field takes a non-zero initial value."
!if(parameters%grid%origin.eq.'nostag') stop 'Trumpet initial data with scalar field not prepared for the nostag grid.'
write(*,*) "Warning: electromagnetic part not evolved (em = 0 hardcoded)."
if(parameters%grid%origin.eq.'nostag') print*, 'Testing nostag scalar field initial data on flat spacetime.'
ai = a; centeri = center; sigmai = sigma
aimi = aim; centerimi = centerim; sigmaimi = sigmaim
aoi = ao; centeroi = centero; sigmaoi = sigmao
!print*, lb, li, ui, ub
dxv=dxv/fac
ui = li + parameters%grid%ncells*fac - 1
if(parameters%grid%origin.eq.'nostag'.or.parameters%grid%origin.eq.'misman') ui = li + parameters%grid%ncells*fac !06-07-2016 misman
ub = ui + parameters%grid%nghost
!print*, lb, li, ui, ub
!!!
!  if (parameters%physics%mass.gt.1e-15) write(*,*) "Mass of the Klein-Gordon field automatically set to zero."

  psii = 1._iwp

  987 continue

  do counter = 1+1+parameters%grid%nghost, size(x)*fac*2-(parameters%grid%nghost*2-1+parameters%grid%nghost), 2
	numr = counter+4
	psii2(counter+1) = (-psii(int((-3 + numr)/2._wp)) + 9*psii(int((-1 + numr)/2._wp)) + 9*psii(int((1 + numr)/2._wp)) - psii(int((3 + numr)/2._wp)))/16._iwp
!print*, psii2(counter+1), psii(int((1 + numr)/2._wp))
  enddo
  do counter = 1,size(x)*fac*2, 2
	numr = counter+4
	psii2(counter+1) = psii(int((numr-1)/2._wp))
  enddo
!print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  do counter = 1,size(x)*fac*2
!print*, psii2(counter), psicompl(counter)
  enddo

!!!!!! Integration of psiA

if(.true.)then !Arr+psiA ! working one
  psiAcompl = 0._iwp
!print*, (1+1+parameters%grid%nghost-2)*2 -1, psiAcompl((1+1+parameters%grid%nghost-2)*2 -1)
!even or odd
  dphikgi = 0
  pikgi = 0
  dimphikgi = 0
  impikgi = 0
  xiv = 0
  prhsi = 0 !crappy fix!! (as many others ...)
  do counter = 1+1+parameters%grid%nghost, int(ubound(psiAcompl,1)/2._iwp)+2 !size(x)*fac
  numr = (counter-2)*2 -1 ! also for nostag
!print*, counter, numr, xiv, psicompl(numr), psiAcompl(numr)
!RK 4th order
  k1i = dxv*prhsi
  xiv = xiv+dxv/2._iwp
  dphikgi = kgf*valuedphiomtbh(xiv,ai,centeri,sigmai)
  pikgi = kgf*valuepiomtbh(xiv,ai,centeri,sigmai,massi,Qi,Kcmci,Ccmci,omegai(numr+1),psicompl(numr+1),parameters%physics%propsign)
  dimphikgi = kgf*valuedphiomtbh(xiv,aimi,centerimi,sigmaimi)
  impikgi = kgf*valuepiomtbh(xiv,aimi,centerimi,sigmaimi,massi,Qi,Kcmci,Ccmci,omegai(numr+1),psicompl(numr+1),parameters%physics%propsign)
  call p2rhsomtpsiA(xiv,psiAcompl(numr-2)+k1i/2._iwp,massi,Qi,Kcmci,Ccmci,qpi,aai,omegai(numr+1),psicompl(numr+1),0._iwp,dphikgi,pikgi,0._iwp,dimphikgi,impikgi,prhsi)
!print*, prhsi
  k2i = dxv*prhsi
  call p2rhsomtpsiA(xiv,psiAcompl(numr-2)+k2i/2._iwp,massi,Qi,Kcmci,Ccmci,qpi,aai,omegai(numr+1),psicompl(numr+1),0._iwp,dphikgi,pikgi,0._iwp,dimphikgi,impikgi,prhsi)
  k3i = dxv*prhsi
  xiv = xiv+dxv/2._iwp
  dphikgi = kgf*valuedphiomtbh(xiv,ai,centeri,sigmai)
  pikgi = kgf*valuepiomtbh(xiv,ai,centeri,sigmai,massi,Qi,Kcmci,Ccmci,omegai(numr+2),psicompl(numr+2),parameters%physics%propsign)
  dimphikgi = kgf*valuedphiomtbh(xiv,aimi,centerimi,sigmaimi)
  impikgi = kgf*valuepiomtbh(xiv,aimi,centerimi,sigmaimi,massi,Qi,Kcmci,Ccmci,omegai(numr+2),psicompl(numr+2),parameters%physics%propsign)
  call p2rhsomtpsiA(xiv,psiAcompl(numr-2)+k3i,massi,Qi,Kcmci,Ccmci,qpi,aai,omegai(numr+2),psicompl(numr+2),0._iwp,dphikgi,pikgi,0._iwp,dimphikgi,impikgi,prhsi)
  k4i = dx*prhsi
  psiAcompl(numr) = psiAcompl(numr-2) + kgf*(k1i+2*k2i+2*k3i+k4i)/6._iwp
  call p2rhsomtpsiA(xiv,psiAcompl(numr),massi,Qi,Kcmci,Ccmci,qpi,aai,omegai(numr+2),psicompl(numr+2),0._iwp,dphikgi,pikgi,0._iwp,dimphikgi,impikgi,prhsi)
!  print*, xiv, psiAcompl(numr), psicompl(numr+2) ! I am not completely happy with outermost values!! there is a negative one!
  enddo
write(*,*)
!odd or even
  dphikgi = 0
  pikgi = 0
  xiv = dxv/2._iwp
  prhsi = 0 !crappy fix!! (as many others ...)
  do counter = 1+1+parameters%grid%nghost, int(ubound(psiAcompl,1)/2._iwp)+2 !size(x)*fac
  numr = (counter-2)*2 ! seems ok, also for nostag
!print*, counter, numr, xiv, psicompl(numr), psiAcompl(numr)
!RK 4th order
  k1i = dxv*prhsi
  xiv = xiv+dxv/2._iwp
  dphikgi = kgf*valuedphiomtbh(xiv,ai,centeri,sigmai)
  pikgi = kgf*valuepiomtbh(xiv,ai,centeri,sigmai,massi,Qi,Kcmci,Ccmci,omegai(numr+1),psicompl(numr+1),parameters%physics%propsign)
  dimphikgi = kgf*valuedphiomtbh(xiv,aimi,centerimi,sigmaimi)
  impikgi = kgf*valuepiomtbh(xiv,aimi,centerimi,sigmaimi,massi,Qi,Kcmci,Ccmci,omegai(numr+1),psicompl(numr+1),parameters%physics%propsign)
  call p2rhsomtpsiA(xiv,psiAcompl(numr-2)+k1i/2._iwp,massi,Qi,Kcmci,Ccmci,qpi,aai,omegai(numr+1),psicompl(numr+1),0._iwp,dphikgi,pikgi,0._iwp,dimphikgi,impikgi,prhsi)
  k2i = dxv*prhsi
  call p2rhsomtpsiA(xiv,psiAcompl(numr-2)+k2i/2._iwp,massi,Qi,Kcmci,Ccmci,qpi,aai,omegai(numr+1),psicompl(numr+1),0._iwp,dphikgi,pikgi,0._iwp,dimphikgi,impikgi,prhsi)
  k3i = dxv*prhsi
  xiv = xiv+dxv/2._iwp
  dphikgi = kgf*valuedphiomtbh(xiv,ai,centeri,sigmai)
  pikgi = kgf*valuepiomtbh(xiv,ai,centeri,sigmai,massi,Qi,Kcmci,Ccmci,omegai(numr+2),psicompl(numr+2),parameters%physics%propsign)
  dimphikgi = kgf*valuedphiomtbh(xiv,aimi,centerimi,sigmaimi)
  impikgi = kgf*valuepiomtbh(xiv,aimi,centerimi,sigmaimi,massi,Qi,Kcmci,Ccmci,omegai(numr+2),psicompl(numr+2),parameters%physics%propsign)
  call p2rhsomtpsiA(xiv,psiAcompl(numr-2)+k3i,massi,Qi,Kcmci,Ccmci,qpi,aai,omegai(numr+2),psicompl(numr+2),0._iwp,dphikgi,pikgi,0._iwp,dimphikgi,impikgi,prhsi)
  k4i = dx*prhsi
  psiAcompl(numr) = psiAcompl(numr-2) + kgf*(k1i+2*k2i+2*k3i+k4i)/6._iwp
  call p2rhsomtpsiA(xiv,psiAcompl(numr),massi,Qi,Kcmci,Ccmci,qpi,aai,omegai(numr+2),psicompl(numr+2),0._iwp,dphikgi,pikgi,0._iwp,dimphikgi,impikgi,prhsi)
!  print*, xiv, psiAcompl(numr), psicompl(numr+2)
  enddo
endif

  psiA(li:ui) = psiAcompl(li-1:ui:2) ! stag ! seems ok for misman too
  psiA(lb:li-1) = psiA(li)
  psiA(ui+1:ub) = psiA(ui)

if (parameters%grid%origin.eq.'nostag') then
  psiA(li:ui) = psiAcompl(li-2:ui:2) !seems ok
  psiA(lb:li-1) = psiA(li)
  psiA(ui+1:ub) = psiA(ui)
endif

!do counter = lbound(psiAcompl,1), ubound(psiAcompl,1)
!print*, psiAcompl(counter)
!enddo
!do counter = lb, ub
!print*, x(counter), psic(counter), psiA(counter)
!enddo

   do counter = 1             , size(x)*fac*2-(parameters%grid%nghost*2-1+parameters%grid%nghost)
write(561,*) dxv/2._iwp*(counter-3), psiAcompl(counter) , psicompl(counter)
   enddo
write(561,*)

!!!!!! Integration of psiAcompl

!06-07-2016 misman: looks ok so far ... 

!if(.false.)then !ooo!

  icond1i = 1.0_iwp !0.9_wp
  icond2i = 1.1_iwp !1.9_wp !+10*a*center*center
if(abs(mass).lt.1d-12) then
  icond1i = 0.9_iwp
  icond2i = 1.08_iwp !1.1_iwp !1.9_wp
endif
  nt = 0
  if(abs(parameters%grid%xmin).gt.1d-12) stop "Leftmost boundary not set at the imposed xmin = 0."
if(parameters%grid%origin.eq.'stag')then
	call epsiomtbh(icond1i,lb,li,ui,ub,0._iwp+dxv/2._iwp,dxv,epsi1i,outi,parameters,omegai,psicompl,psiAcompl)
	call epsiomtbh(icond2i,lb,li,ui,ub,0._iwp+dxv/2._iwp,dxv,epsi2i,outi,parameters,omegai,psicompl,psiAcompl)
elseif(parameters%grid%origin.eq.'nostag')then
	call epsiomtnostagbh(icond1i,lb,li,ui,ub,0._iwp+dxv,dxv,epsi1i,outi,parameters,omegai,psicompl,psiAcompl)
	call epsiomtnostagbh(icond2i,lb,li,ui,ub,0._iwp+dxv,dxv,epsi2i,outi,parameters,omegai,psicompl,psiAcompl)
elseif(parameters%grid%origin.eq.'misman')then
	call epsiomtmismanbh(icond1i,lb,li,ui,ub,0._iwp+dxv/2._iwp,dxv,epsi1i,outi,parameters,omegai,psicompl,psiAcompl)
	call epsiomtmismanbh(icond2i,lb,li,ui,ub,0._iwp+dxv/2._iwp,dxv,epsi2i,outi,parameters,omegai,psicompl,psiAcompl)
endif

  658 continue
  nt = nt + 1
  if(mod(nt,50).eq.0) print*, "Current error function at scri = ", epsi3i
  if(nt.gt.200) stop "nt larger than 200!"
!  print*, epsi1, epsi2
! write(555,*) icond1, epsi1
! write(555,*) icond2, epsi2
  if (epsi1i*epsi2i.gt.0) then
  write(*,*) "Initial values for the bisection for psi not correct: "
  write(*,'(2(f10.4,a,e10.4,a))') icond1i, ' -> ', epsi1i, ',  ', icond2i, ' -> ', epsi2i
  stop
  endif
  icond3i = (icond2i+icond1i)/2._iwp
if(parameters%grid%origin.eq.'stag')then
	call epsiomtbh(icond3i,lb,li,ui,ub,0._iwp+dxv/2._iwp,dxv,epsi3i,outi,parameters,omegai,psicompl,psiAcompl)
elseif(parameters%grid%origin.eq.'nostag')then
	call epsiomtnostagbh(icond3i,lb,li,ui,ub,0._iwp+dxv,dxv,epsi3i,outi,parameters,omegai,psicompl,psiAcompl)
elseif(parameters%grid%origin.eq.'misman')then
	call epsiomtmismanbh(icond3i,lb,li,ui,ub,0._iwp+dxv/2._iwp,dxv,epsi3i,outi,parameters,omegai,psicompl,psiAcompl)
endif
!  print*, nt!
!stop
  if (abs(epsi3i).gt.1d-10) then !8 !12
    if (epsi1i*epsi3i.gt.0) then
 icond1i = icond3i
 epsi1i = epsi3i
    else
 icond2i = icond3i
 epsi2i = epsi3i
    endif
    goto 658
  endif
  write(*,*) '-- number of iterations = ', nt, epsi3i, psi(li), psi(ui)
  write(*,*) 'selected initial condition = ', icond3i
 !close(555)
!!!
dxv=dxv*fac
ui = li + parameters%grid%ncells - 1
if(parameters%grid%origin.eq.'nostag'.or.parameters%grid%origin.eq.'misman') ui = li + parameters%grid%ncells
ub = ui + parameters%grid%nghost
!print*, lb, li, ui, ub
psi(li:ui) = outi(li+int((fac-1)/2._wp):ui:fac)
psii(li:ui) = outi(li+int((fac-1)/2._wp):ui:fac)
!if(parameters%grid%origin.eq.'nostag') psi(li:ui) = out(li:ui:fac)
!!!

do counter = lb, ub
	write(560,*) x(counter), psi(counter), out(counter)
enddo

 if (parameters%grid%origin.eq.'nostag') then 
	select case (parameters%moldef%deriv_method)
	   case ('c2')
  		psi(li) = (4.0_wp*psi(li+1) - psi(li+2))/3.0_wp
	   case ('c4')
  		psi(li) = (48._wp*psi(li+1) - 36._wp*psi(li+2) + 16._wp*psi(li+3) - 3._wp*psi(li+4))/25.0_wp
	   case ('c6')
  		psi(li) = (360._wp*psi(li+1) - 450._wp*psi(li+2) + 400._wp*psi(li+3) - 225._wp*psi(li+4) + 72._wp*psi(li+5) - 10._wp*psi(li+6))/147.0_wp
	   case ('c8')
  		psi(li) = (6720._wp*psi(li+1) - 11760._wp*psi(li+2) + 15680._wp*psi(li+3) - 14700._wp*psi(li+4) + 9408._wp*psi(li+5) - 3920._wp*psi(li+6) + 960._wp*psi(li+7) - 105._wp*psi(li+8))/2283.0_wp
	end select
	do counter = 1, abs(parameters%grid%nghost)
  		psi(lb + counter - 1) = + psi(li + abs(parameters%grid%nghost) - counter + 1)
		psii(lb + counter - 1) = + psii(li + abs(parameters%grid%nghost) - counter + 1)
		psi(ui+counter) = 5*psi(ui+counter-1)-10*psi(ui+counter-2)+10*psi(ui+counter-3)-5*psi(ui+counter-4)+psi(ui+counter-5)
		psii(ui+counter) = 5*psii(ui+counter-1)-10*psii(ui+counter-2)+10*psii(ui+counter-3)-5*psii(ui+counter-4)+psii(ui+counter-5)
	enddo

    elseif (parameters%grid%origin.eq.'stag') then
	do counter = 1, abs(parameters%grid%nghost)
  		psi(lb + counter - 1) = + psi(li + abs(parameters%grid%nghost) - counter)
		psii(lb + counter - 1) = + psii(li + abs(parameters%grid%nghost) - counter)
		psi(ui+counter) = 5*psi(ui+counter-1)-10*psi(ui+counter-2)+10*psi(ui+counter-3)-5*psi(ui+counter-4)+psi(ui+counter-5)
		psii(ui+counter) = 5*psii(ui+counter-1)-10*psii(ui+counter-2)+10*psii(ui+counter-3)-5*psii(ui+counter-4)+psii(ui+counter-5)
	enddo

	elseif (parameters%grid%origin.eq.'misman') then !06-07-2016 misman ! looks the same as the stag case ... 
	do counter = 1, abs(parameters%grid%nghost)
  		psi(lb + counter - 1) = + psi(li + abs(parameters%grid%nghost) - counter)
		psii(lb + counter - 1) = + psii(li + abs(parameters%grid%nghost) - counter)
		psi(ui+counter) = 5*psi(ui+counter-1)-10*psi(ui+counter-2)+10*psi(ui+counter-3)-5*psi(ui+counter-4)+psi(ui+counter-5)
		psii(ui+counter) = 5*psii(ui+counter-1)-10*psii(ui+counter-2)+10*psii(ui+counter-3)-5*psii(ui+counter-4)+psii(ui+counter-5)
	enddo

    else
	stop 'Is the grid staggered or not?' 
	!continue
    end if

do counter = lb, ub
!working here print*, counter, x(counter), psi(counter)
	write(559,*) x(counter), psi(counter)
enddo
write(559,*)

!endif !oooo!
!psi=1 !oooo!

  u(1,:) = u(1,:)/psi**4

  u(3,:) = (u(3,:)+psiA)/psi**6

  u(9,:) = noise*u(9,:) + kgf*(a*exp(-(x**2-center**2)**2/(4._wp*sigma**4)))/omega ! moved from above for rescaled phikg -- now rephi
  
  u(10,:) = -(-a*exp(-(x**2-center**2)**2/(4._wp*sigma**4))*x*(x**2-center**2)/sigma**4) * ( Ccmc*psic**3/x**2 + Kcmc*x/3._iwp )*(-propsign/psi**6-1)/omega ! for rescaled pikg (pphi=0) -- now repi
  !to change direction of movement of phi, change sign in front of 1/psi**6 and in valuepiomtbh!!!
  u(11,:) = noise*u(11,:) + kgf*(aim*exp(-(x**2-centerim**2)**2/(4._wp*sigmaim**4)))/omega ! imaginary parts of scalar field
  u(12,:) = -(-aim*exp(-(x**2-centerim**2)**2/(4._wp*sigmaim**4))*x*(x**2-centerim**2)/sigmaim**4) * ( Ccmc*psic**3/x**2 + Kcmc*x/3._iwp )*(-propsign/psi**6-1)/omega ! imaginary parts of scalar field

endif ! solving for the perturbation (not needed for Cowling)

  u(14,:) = Q*u(1,:)**(1.5_wp)/x**2 ! Er

 if (parameters%grid%origin.eq.'nostag') then !setting to zero initial value of the scalar field variables at scri: necessary when the rescaled variables are used
 print*, 'Warning: setting central value of pikg manually.' ! in the following line
 u(10,li) = -(-a*exp(-(x(li)**2-center**2)**2/(4._wp*sigma**4))*x(li)*(x(li)**2-center**2)/sigma**4) * ( 0 + Kcmc*x(li)/3._iwp )*(-propsign/psi(li)**6-1)/omega(li) ! for rescaled pikg (pphi=0)
 u(9,ui) = 0
 u(10,ui) = 0
 u(11,ui) = 0
 u(12,ui) = 0
 endif
 
if (parameters%grid%origin.eq.'misman') then !setting to zero initial value of the scalar field variables at scri: necessary when the rescaled variables are used !06-07-2016 misman
 u(9,ui) = 0
 u(10,ui) = 0
 u(11,ui) = 0
 u(12,ui) = 0
 endif


print*, 'inner chi = ', u(1,li), ',   outer chi = ', u(1,ui)
print*, 'inner Arr = ', u(3,li), ',   outer Arr = ', u(3,ui)

!linear! endif
!! Initial data KGF hypcomp with BH

    endif

include '../source.f90' ! to calculate the constraints

  case ('only_gauge')

!  if (.false.) then !attempt evolving only the gauge conditions, trying to determine their stationary solution
  u = 0
do i = lbound(x,1), ubound(x,1)
beta(i)=psic(i)
enddo
  if (parameters%moldef%deriv_method == 'c8') then
     call vdiff_c8(beta, dpsic, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
     call vdiff_c6(beta, dpsic, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call vdiff_c4(beta, dpsic, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
     call vdiff(beta, dpsic, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if

!  call vdiff_c4(psic, dpsic, x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
do i = lbound(x,1), ubound(x,1)
  u(3,i) = dpsic(i)/omega(i)-psic(i)*domega(i)/omega(i)**2
enddo
  u(7,:) = 0*u(7,:) + Ccmc*psic**3/x**2 + Kcmc*x/3._wp ! this was giving me so much trouble ... without the first term it would return just 0 ... why? no idea ...
  do counter = 1, abs(parameters%grid%nghost)
  		u(3,lbound(x,1) + counter - 1) = u(3,lbound(x,1) + 2*abs(parameters%grid%nghost) - counter)
  		u(3,ubound(x,1) - counter + 1) = u(3,ubound(x,1) - 2*abs(parameters%grid%nghost) + counter)
  enddo
!  do i = lbound(x,1), ubound(x,1)
 ! write(888,*) x(i), dpsic(i), u(7,i), u(3,i)
  !enddo
  
!  endif

include '../source.f90' ! to calculate the constraints  

  case ('perturb_relaxed_trumpet')

stop "Not updated to include perturbations other than by a real scalar field."

open(13, FILE=trim(file_prefix)//'last_data.txt', STATUS='OLD')
read(13,*) n, parameters%timer%t
parameters%timer%end_time = parameters%timer%end_time + parameters%timer%t
!if(n.ne.(2*(size(u(1,:))-2*parameters%grid%nghost)+2*parameters%grid%nghost)) then 
if(n.ne.size(u(1,:))) then 
write(*,*) n, size(u(1,:))
stop "wrong number of data points from relaxed trumpet"
end if
!write(*,*) (parameters%grid%ncells+1)-n,size(u(1,:))-n
write(*,*) lbound(u,2), ubound(u,2)
do i = lbound(u,2), ubound(u,2)
read(13,'(16(e35.24e3))') x(i), (u(j,i),j=lbound(u,1), ubound(u,1)) !! added xdo(i), 
!write(*,*) xdo(i), (udo(j,i),j=lbound(udo,1), ubound(udo,1)) !looks alright
enddo
 close(13)

!print*, xdo(lbound(xdo)),xdo(ubound(xdo)), xdo(lbound(udo,2)),xdo(ubound(udo,2)) !looks alright

  icond1 = 0.9_wp 
  icond2 = 1.3_wp 
if(abs(mass).lt.1d-12) then
  icond1 = 0.9_wp
  icond2 = 1.9_wp
endif
  nt = 0
  if(abs(parameters%grid%xmin).gt.1d-12) stop "Leftmost boundary not set at the imposed xmin = 0."
if(parameters%grid%origin.eq.'stag')then
	call epsirelaxtrumpet(icond1,a,center,sigma,lb,li,ui,ub,x(lbound(u,2)),x(ubound(u,2)),0._wp+dx/2._wp,dx,epsi1,out,parameters,x,u)
			print*, icond1, epsi1
	call epsirelaxtrumpet(icond2,a,center,sigma,lb,li,ui,ub,x(lbound(u,2)),x(ubound(u,2)),0._wp+dx/2._wp,dx,epsi2,out,parameters,x,u)
			print*, icond2, epsi2
else
 stop "Include other grid options!"
endif

  668 continue
  nt = nt + 1
  if(mod(nt,50).eq.0) print*, "Current error function at scri = ", epsi3
  if(nt.gt.200) stop "nt larger than 200!"
!  print*, epsi1, epsi2, icond1, icond2
! write(555,*) icond1, epsi1
! write(555,*) icond2, epsi2
  if (epsi1*epsi2.gt.0) then
  write(*,*) "Initial values for the bisection for psi (relaxed trumpet) not correct: "
  write(*,'(4(f10.4,a))') icond1, ' -> ', epsi1, ',  ', icond2, ' -> ', epsi2
  stop
  endif
  icond3 = (icond2+icond1)/2._wp
if(parameters%grid%origin.eq.'stag')then
	call epsirelaxtrumpet(icond3,a,center,sigma,lb,li,ui,ub,x(lbound(u,2)),x(ubound(u,2)),0._wp+dx/2._wp,dx,epsi3,out,parameters,x,u)
	print*, icond3, epsi3
endif
!  print*, nt!
!stop
  if (abs(epsi3).gt.1d-10) then !8 !12
    if (epsi1*epsi3.gt.0) then
 icond1 = icond3
 epsi1 = epsi3
    else
 icond2 = icond3
 epsi2 = epsi3
    endif
    goto 668
  endif
  write(*,*) '-- number of iterations = ', nt, epsi3, out(li), out(ui)
  write(*,*) 'selected initial condition = ', icond3
 !close(555)
!!!
! now fac is 1 anyway dxv=dxv*fac
ui = li + parameters%grid%ncells - 1
if(parameters%grid%origin.eq.'nostag'.or.parameters%grid%origin.eq.'misman') ui = li + parameters%grid%ncells
ub = ui + parameters%grid%nghost
!print*, lb, li, ui, ub
psi(li:ui) = outi(li+int((fac-1)/2._wp):ui:fac)
!!!

do counter = lb, ub
	write(560,*) x(counter), psi(counter), out(counter)
enddo

 if (parameters%grid%origin.eq.'nostag') then 
	select case (parameters%moldef%deriv_method)
	   case ('c2')
  		psi(li) = (4.0_wp*psi(li+1) - psi(li+2))/3.0_wp
	   case ('c4')
  		psi(li) = (48._wp*psi(li+1) - 36._wp*psi(li+2) + 16._wp*psi(li+3) - 3._wp*psi(li+4))/25.0_wp
	   case ('c6')
  		psi(li) = (360._wp*psi(li+1) - 450._wp*psi(li+2) + 400._wp*psi(li+3) - 225._wp*psi(li+4) + 72._wp*psi(li+5) - 10._wp*psi(li+6))/147.0_wp
	   case ('c8')
  		psi(li) = (6720._wp*psi(li+1) - 11760._wp*psi(li+2) + 15680._wp*psi(li+3) - 14700._wp*psi(li+4) + 9408._wp*psi(li+5) - 3920._wp*psi(li+6) + 960._wp*psi(li+7) - 105._wp*psi(li+8))/2283.0_wp
	end select
	do counter = 1, abs(parameters%grid%nghost)
  		psi(lb + counter - 1) = + psi(li + abs(parameters%grid%nghost) - counter + 1)
		psi(ui+counter) = 5*psi(ui+counter-1)-10*psi(ui+counter-2)+10*psi(ui+counter-3)-5*psi(ui+counter-4)+psi(ui+counter-5)
	enddo

    elseif (parameters%grid%origin.eq.'stag') then
	do counter = 1, abs(parameters%grid%nghost)
  		psi(lb + counter - 1) = + psi(li + abs(parameters%grid%nghost) - counter)
		psi(ui+counter) = 5*psi(ui+counter-1)-10*psi(ui+counter-2)+10*psi(ui+counter-3)-5*psi(ui+counter-4)+psi(ui+counter-5)
	enddo

	elseif (parameters%grid%origin.eq.'misman') then !06-07-2016 misman ! looks the same as the stag case ... 
	do counter = 1, abs(parameters%grid%nghost)
  		psi(lb + counter - 1) = + psi(li + abs(parameters%grid%nghost) - counter)
		psi(ui+counter) = 5*psi(ui+counter-1)-10*psi(ui+counter-2)+10*psi(ui+counter-3)-5*psi(ui+counter-4)+psi(ui+counter-5)
	enddo

    else
	stop 'Is the grid staggered or not?' 
	!continue
    end if

do counter = lb, ub
!working here print*, counter, x(counter), psi(counter)
	write(559,*) x(counter), psi(counter)
enddo
write(559,*)

  u = udo(:,::2) !! check if indices correct

  u(1,:) = u(1,:)/psi**4

!! only considering psi now  

print*, 'inner chi = ', u(1,li), ',   outer chi = ', u(1,ui)
print*, 'inner Arr = ', u(3,li), ',   outer Arr = ', u(3,ui)


include '../source.f90' ! to calculate the constraints
 


   case default

stop 'No initial data specified.'

   end select

end subroutine initialdata






subroutine setboundaries(vect, x, parameters)
  use numservicef90
  implicit none

    real(kind=wp), dimension(:,:),intent(inout) :: vect
    real(kind=wp), dimension(:),intent(inout) :: x
    type(parameters_type), intent(in)    :: parameters
!*
  real(kind=wp), dimension(size(vect,1), size(vect,2))   :: u, gauge !,intent(inout)
  real(kind=wp), dimension(size(u,1), size(u,2))  :: constr !,intent(out)
!  real(kind=wp), dimension(:),intent(in)    :: x
!  type(parameters_type), intent(inout) :: parameters
  real(kind=wp), dimension(size(u,1), size(u,2)):: rhs !, intent(out)
!
  real(kind=wp), dimension(size(u,1), size(u,2)) :: du, ddu, u_diss, dru, druin
  real(kind=wp), dimension(size(x)) :: chi, grr, gtt, Arr, K, Lambdar, alpha, betar, Br, derchi, dergrr, dergtt, deralpha, derbetar, phikg, pikg, Theta, Zr, Lambdazr, RR, TT, rephi, repi, imphi, impi, mEr, mAr, mPhi, mPsi
  real(kind=wp), dimension(size(x)) :: dchi, dgrr, dgtt, dArr, dK, dLambdar, dalpha, dbetar, dBr, dderchi, ddergrr, ddergtt, dderalpha, dderbetar, dphikg, dpikg, dTheta, dZr, dLambdazr, dRR, dTT, drephi, drepi, dimphi, dimpi, dmEr, dmAr, dmPhi, dmPsi
  real(kind=wp), dimension(size(x)) :: drchi, drgrr, drgtt, drArr, drK, drLambdar, dralpha, drbetar, drBr, drphikg, drpikg, drTheta, drZr, drRR, drTT, drrephi, drrepi, drimphi, drimpi, drmEr, drmAr, drmPhi, drmPsi ! for advection terms
  real(kind=wp), dimension(size(x)) :: ddchi, ddgrr, ddgtt, ddArr, ddK, ddLambdar, ddalpha, ddbetar, ddBr, ddphikg, ddpikg, ddTheta, ddZr, ddrephi, ddrepi, ddimphi, ddimpi, ddmEr, ddmAr, ddmPhi, ddmPsi
!  real(kind=wp), dimension(size(x)) :: auxmass, auxp, auxpp, div, squ, sqar, sqar2, sqar3, gsqarba, gsqarba2, ugsqarb, usqarb, gsqarb
  real(kind=wp) :: v, eta, lambda, mu, eps, mass, cosm, cosmk, aa, cuplog, charm, xmax, xscri, coealphah, coealphal, coebetar, coeLambdar, coeZr, pphi, trafK, Kcmc, Ccmc, k1, xexc
  integer	:: evollapse, evolshift, last, nro, nr0, kgf, null, numr, z4, z4c, omt, em, gr !counter,

  real(kind=wp), dimension(size(x)) :: auxmass, auxp, auxpp, div, squ, omega, domega, ddomega, dddomega, root, root2, coefom, ders, der ! sqar, sqar2, sqar3, gsqarba, gsqarba2, ugsqarb, usqarb, gsqarb,
  real(kind=wp), dimension(size(x)) :: dersdalpha, dersdgtt, dersdchi, dersArr,dersK, dersLambdar, dersalphau, dedersgttchi, ddersgttchi
  real(kind=wp), dimension(size(x)) :: grrb, dgrrb, ddgrrb, gttb, dgttb, ddgttb, coefalpha, dcoefalpha, ddcoefalpha, coefgrr, dcoefgrr, ddcoefgrr, coefgrrb, dcoefgrrb, ddcoefgrrb, nuplog, nconst
  real(kind=wp) :: phiscri, piscri, alphascri, grrscri, chiscri, betarscri
!*
    integer :: lb, ub, li, ui, counter
    real(kind=wp) :: dx, coeff
  real(kind=wp), dimension(size(x)) :: transition

  integer, parameter :: n2v = 18
!*
  cuplog = parameters%slice%cuplog
  charm = parameters%slice%charm
  evollapse = parameters%slice%evollapse
  evolshift = parameters%slice%evolshift
  v = parameters%physics%v
  eta = parameters%physics%eta
  lambda = parameters%physics%lambda
  mu = parameters%physics%mu
  aa = parameters%physics%aa
  mass = parameters%physics%mass
  cosm = parameters%physics%cosm
  gr = parameters%physics%gr
  kgf = parameters%physics%kgf
  em = parameters%physics%em
  omt = parameters%physics%omt
  cosmk = cosm/3._wp
  xmax = parameters%grid%xmax
  xscri = parameters%physics%xscri
  xexc = parameters%physics%xexc
  k1 = parameters%physics%k1
  pphi = parameters%physics%pphi
  Kcmc = parameters%physics%Kcmc
  coealphah = parameters%physics%coealphah

  if (aa.lt.0) aa = -3/Kcmc ! not really needed - have to clean up

  eps = parameters%moldef%dissipation_eps
!*

!*

  dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / dble(parameters%grid%ncells)
!  if (parameters%grid%origin.eq.'stag') dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / (dble(parameters%grid%ncells+1)) !change for staggered grid ! RUINS CONVERGENCE!!!!!!
  if (parameters%grid%origin.eq.'misman') dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / (dble(parameters%grid%ncells)+0.5_wp) !06-07-2016 misman

!*

    ! executable statements

    lb = lbound(vect, 2)
    ub = ubound(vect, 2)

    li = lb + parameters%grid%nghost
    ui = ub - parameters%grid%nghost

!    dx = abs(parameters%grid%xmax - parameters%grid%xmin) / dble(parameters%grid%ncells)

    select case (parameters%moldef%boundary_type)

!bc
    case ('zero')
  ! From li to lb and from ui to ub because
  ! ui + li - lb = ub - nghost + lb + nghost - lb = ub
  ! That sets also the ghost points to zero.
  do counter = 0, li - lb
vect(:, li-counter) = 0.0_wp
vect(:, ui+counter) = 0.0_wp
  end do

!bc
    case ('periodic')
  do counter = 0, abs(parameters%grid%nghost)
vect(:, ui + counter) = vect(:, li + counter)
vect(:, li - counter) = vect(:, ui - counter)
  enddo

!bc
    case ('reflecting')
  vect(1, lb) = 0.0_wp
  vect(1, ub) = 0.0_wp

  vect(2, lb) = vect(2, li)
  vect(2, ub) = vect(2, ui)

  vect(3, lb) = 0.0_wp
  vect(3, ub) = 0.0_wp

!bc
    case('sommerfeld')
  vect(2, li) =   vect(1, li)
  vect(2, ui) = - vect(1, ui)

!bc
    case('compact')
  ! Inner boundary at the center
  vect(2, li) = 0d0

!bc
    case('spherical_hyperal')
  ! Inner boundary at the center
  vect(2, li) = 0d0
  ! Outgoing at the outer boundary
  vect(2, ui) = - vect(1, ui) - vect(3,ui)!/parameters%grid%xmax

!bc
    case('spherical')
  ! Inner boundary at the center
  vect(2, li) = 0d0

  ! Set the ingoing characteristic to zero at the outer boundary
  vect(2, ui) = - vect(1, ui) - vect(3,ui)/parameters%grid%xmax  !psi = -pi - phi

!bc
    case('4waver')
	do counter = 1, abs(parameters%grid%nghost)
  		vect(1, lb + counter - 1) = + vect(1, li + abs(parameters%grid%nghost) - counter) ! +
  	    	!vect(1, ub - counter + 1) = 0.0_wp
	   	!vect(2, lb + counter - 1) = - vect(2, li + abs(parameters%grid%nghost) - counter) ! -
   	 	!vect(2, ub - counter + 1) = 0.0_wp
  		vect(2, lb + counter - 1) = + vect(2, li + abs(parameters%grid%nghost) - counter) ! +
	 	!vect(3, ub - counter + 1) = 0.0_wp
		!vect(1,ub - counter + 1) = 0.0_wp
		!vect(2,ub - counter + 1) = vect(1,ub - counter + 1)
	enddo
	!vect(2, ui) = - vect(1, ui) - vect(3,ui)!/parameters%grid%xmax

!bc
    case('gbssn')

    !Non-staggered grid
    if (parameters%grid%origin.eq.'nostag') then

! if(.true.) then
	! at xmin
	select case (parameters%moldef%deriv_method)
	   case ('c2')
  		vect(:, li) = (4.0_wp*vect(:, li+1) - vect(:, li+2))/3.0_wp
  		vect(5, li) = 0.0_wp
  		vect(7, li) = 0.0_wp
  		vect(13:14, li) = 0.0_wp
  		!vect(n2v+1:n2v+4, li) = 0.0_wp
	   case ('c4')
  		vect(:, li) = (48._wp*vect(:,li+1) - 36._wp*vect(:,li+2) + 16._wp*vect(:,li+3) - 3._wp*vect(:,li+4))/25.0_wp
  		vect(5, li) = 0.0_wp
  		vect(7, li) = 0.0_wp
  		vect(13:14, li) = 0.0_wp
  		!vect(n2v+1:n2v+4, li) = 0.0_wp
	   case ('c6')
  		vect(:, li) = (360._wp*vect(:,li+1) - 450._wp*vect(:,li+2) + 400._wp*vect(:,li+3) - 225._wp*vect(:,li+4) + 72._wp*vect(:,li+5) - 10._wp*vect(:,li+6))/147.0_wp
  		vect(5, li) = 0.0_wp
  		vect(7, li) = 0.0_wp
  		vect(13:14, li) = 0.0_wp
  		!vect(n2v+1:n2v+4, li) = 0.0_wp
	   case ('c8')
  		vect(:, li) = (6720._wp*vect(:,li+1) - 11760._wp*vect(:,li+2) + 15680._wp*vect(:,li+3) - 14700._wp*vect(:,li+4) + 9408._wp*vect(:,li+5) - 3920._wp*vect(:,li+6) + 960._wp*vect(:,li+7) - 105._wp*vect(:,li+8))/2283.0_wp
  		vect(5, li) = 0.0_wp
  		vect(7, li) = 0.0_wp
  		vect(13:14, li) = 0.0_wp
  		!vect(n2v+1:n2v+4, li) = 0.0_wp
	end select
	! at xmax - wrong!!!
!stop 'Wrong conditions at the outer boundary!'
!endif
	do counter = 1, abs(parameters%grid%nghost)
  		vect(:, lb + counter - 1) = + vect(:, li + abs(parameters%grid%nghost) - counter + 1)
  		vect(5, lb + counter - 1) = - vect(5, li + abs(parameters%grid%nghost) - counter + 1)
  		vect(7, lb + counter - 1) = - vect(7, li + abs(parameters%grid%nghost) - counter + 1)
  		vect(13:14, lb + counter - 1) = - vect(13:14, li + abs(parameters%grid%nghost) - counter + 1)
  		!vect(n2v+1:n2v+4, lb + counter - 1) = - vect(n2v+1:n2v+4, li + abs(parameters%grid%nghost) - counter + 1)
	enddo


if(.true.)then ! copied from the stag case - update!!!
! overflow boundary of 0509119
	! at xmax
!	select case (parameters%moldef%deriv_method)
	select case (parameters%grid%nghost)
! the variable's indices changed, so those derived twice need to be relabeled -- for now using the higher order outflow condition for all variables 2024/02/23
!	   case ('c2')
	   case(1)
		!Variables derived only once - (hD_-)^2
  		vect(:, ub-0) = 2*vect(:, ui-0)-vect(:, ui-1)
 		vect(:, ub-0) = 3*vect(:, ui-0)-3*vect(:, ui-1)+vect(:, ui-2) ! apparently very similar
		!Variables derived twice - (hD_-)^3
!  		vect(1:3, ub-0) = 3*vect(1:3, ui-0)-3*vect(1:3, ui-1)+vect(1:3, ui-2)
!  		vect(7:8, ub-0) = 3*vect(7:8, ui-0)-3*vect(7:8, ui-1)+vect(7:8, ui-2)
!  		vect(10, ub-0) = 3*vect(10, ui-0)-3*vect(10, ui-1)+vect(10, ui-2)
!	   case ('c4')
	   case(2)
		!Variables derived only once - (hD_-)^4

!  		vect(:, ub-3) = 4*vect(:, ui-2)-6*vect(:, ui-3)+4*vect(:, ui-4)-vect(:, ui-5)
 ! 		vect(:, ub-2) = 4*vect(:, ui-1)-6*vect(:, ui-2)+4*vect(:, ui-3)-vect(:, ui-4)

  		vect(:, ub-1) = 4*vect(:, ui-0)-6*vect(:, ui-1)+4*vect(:, ui-2)-vect(:, ui-3)
  		vect(:, ub-0) = 4*vect(:, ui+1)-6*vect(:, ui-0)+4*vect(:, ui-1)-vect(:, ui-2)
 		vect(:, ub-1) = 5*vect(:, ui-0)-10*vect(:, ui-1)+10*vect(:, ui-2)-5*vect(:, ui-3)+vect(:, ui-4)
 		vect(:, ub-0) = 5*vect(:, ui+1)-10*vect(:, ui-0)+10*vect(:, ui-1)-5*vect(:, ui-2)+vect(:, ui-3)
		!Variables derived twice - (hD_-)^5

  !		vect(1:3, ub-3) = 5*vect(1:3, ui-2)-10*vect(1:3, ui-3)+10*vect(1:3, ui-4)-5*vect(1:3, ui-5)+vect(1:3, ui-6)
   !    		vect(1:3, ub-2) = 5*vect(1:3, ui-1)-10*vect(1:3, ui-2)+10*vect(1:3, ui-3)-5*vect(1:3, ui-4)+vect(1:3, ui-5)
    !   		vect(7:8, ub-3) = 5*vect(7:8, ui-2)-10*vect(7:8, ui-3)+10*vect(7:8, ui-4)-5*vect(7:8, ui-5)+vect(7:8, ui-6)
!  		vect(7:8, ub-2) = 5*vect(7:8, ui-1)-10*vect(7:8, ui-2)+10*vect(7:8, ui-3)-5*vect(7:8, ui-4)+vect(7:8, ui-5)

!  		vect(1:3, ub-1) = 5*vect(1:3, ui-0)-10*vect(1:3, ui-1)+10*vect(1:3, ui-2)-5*vect(1:3, ui-3)+vect(1:3, ui-4)
!  		vect(1:3, ub-0) = 5*vect(1:3, ui+1)-10*vect(1:3, ui-0)+10*vect(1:3, ui-1)-5*vect(1:3, ui-2)+vect(1:3, ui-3)
!  		vect(7:8, ub-1) = 5*vect(7:8, ui-0)-10*vect(7:8, ui-1)+10*vect(7:8, ui-2)-5*vect(7:8, ui-3)+vect(7:8, ui-4)
!  		vect(7:8, ub-0) = 5*vect(7:8, ui+1)-10*vect(7:8, ui-0)+10*vect(7:8, ui-1)-5*vect(7:8, ui-2)+vect(7:8, ui-3)
!  		vect(10, ub-1) = 5*vect(10, ui-0)-10*vect(10, ui-1)+10*vect(10, ui-2)-5*vect(10, ui-3)+vect(10, ui-4)
!  		vect(10, ub-0) = 5*vect(10, ui+1)-10*vect(10, ui-0)+10*vect(10, ui-1)-5*vect(10, ui-2)+vect(10, ui-3)
!	   case ('c6')
	   case(3)
!stop 'No conditions at the outer boundary!'
		!Variables derived only once - (hD_-)^6
  		vect(:, ub-2) = - vect(:,-5 + ui) + 6*vect(:,-4 + ui) - 15*vect(:,-3 + ui) + 20*vect(:,-2 + ui) - 15*vect(:,-1 + ui) + 6*vect(:,ui)
  		vect(:, ub-1) = - vect(:,-4 + ui) + 6*vect(:,-3 + ui) - 15*vect(:,-2 + ui) + 20*vect(:,-1 + ui) - 15*vect(:,ui) + 6*vect(:,1 + ui)
  		vect(:, ub-0) = - vect(:,-3 + ui) + 6*vect(:,-2 + ui) - 15*vect(:,-1 + ui) + 20*vect(:,ui) - 15*vect(:,1 + ui) + 6*vect(:,2 + ui)
  		vect(:, ub-2) = vect(:,-6 + ui) - 7*vect(:,-5 + ui) + 21*vect(:,-4 + ui) - 35*vect(:,-3 + ui) + 35*vect(:,-2 + ui) - 21*vect(:,-1 + ui) + 7*vect(:,ui)
  		vect(:, ub-1) = vect(:,-5 + ui) - 7*vect(:,-4 + ui) + 21*vect(:,-3 + ui) - 35*vect(:,-2 + ui) + 35*vect(:,-1 + ui) - 21*vect(:,ui) + 7*vect(:,1 + ui)
  		vect(:, ub-0) = vect(:,-4 + ui) - 7*vect(:,-3 + ui) + 21*vect(:,-2 + ui) - 35*vect(:,-1 + ui) + 35*vect(:,ui) - 21*vect(:,1 + ui) + 7*vect(:,2 + ui)
		!Variables derived twice - (hD_-)^7
!		vect(1:3, ub-2) = vect(1:3,-6 + ui) - 7*vect(1:3,-5 + ui) + 21*vect(1:3,-4 + ui) - 35*vect(1:3,-3 + ui) + 35*vect(1:3,-2 + ui) - 21*vect(1:3,-1 + ui) + 7*vect(1:3,ui)
!  		vect(1:3, ub-1) = vect(1:3,-5 + ui) - 7*vect(1:3,-4 + ui) + 21*vect(1:3,-3 + ui) - 35*vect(1:3,-2 + ui) + 35*vect(1:3,-1 + ui) - 21*vect(1:3,ui) + 7*vect(1:3,1 + ui)
!  		vect(1:3, ub-0) = vect(1:3,-4 + ui) - 7*vect(1:3,-3 + ui) + 21*vect(1:3,-2 + ui) - 35*vect(1:3,-1 + ui) + 35*vect(1:3,ui) - 21*vect(1:3,1 + ui) + 7*vect(1:3,2 + ui)
!		vect(7:8, ub-2) = vect(7:8,-6 + ui) - 7*vect(7:8,-5 + ui) + 21*vect(7:8,-4 + ui) - 35*vect(7:8,-3 + ui) + 35*vect(7:8,-2 + ui) - 21*vect(7:8,-1 + ui) + 7*vect(7:8,ui)
!  		vect(7:8, ub-1) = vect(7:8,-5 + ui) - 7*vect(7:8,-4 + ui) + 21*vect(7:8,-3 + ui) - 35*vect(7:8,-2 + ui) + 35*vect(7:8,-1 + ui) - 21*vect(7:8,ui) + 7*vect(7:8,1 + ui)
!  		vect(7:8, ub-0) = vect(7:8,-4 + ui) - 7*vect(7:8,-3 + ui) + 21*vect(7:8,-2 + ui) - 35*vect(7:8,-1 + ui) + 35*vect(7:8,ui) - 21*vect(7:8,1 + ui) + 7*vect(7:8,2 + ui)
!		vect(10, ub-2) = vect(10,-6 + ui) - 7*vect(10,-5 + ui) + 21*vect(10,-4 + ui) - 35*vect(10,-3 + ui) + 35*vect(10,-2 + ui) - 21*vect(10,-1 + ui) + 7*vect(10,ui)
!  		vect(10, ub-1) = vect(10,-5 + ui) - 7*vect(10,-4 + ui) + 21*vect(10,-3 + ui) - 35*vect(10,-2 + ui) + 35*vect(10,-1 + ui) - 21*vect(10,ui) + 7*vect(10,1 + ui)
!  		vect(10, ub-0) = vect(10,-4 + ui) - 7*vect(10,-3 + ui) + 21*vect(10,-2 + ui) - 35*vect(10,-1 + ui) + 35*vect(10,ui) - 21*vect(10,1 + ui) + 7*vect(10,2 + ui)

!	   case ('c8')
	   case(4)
!stop 'No conditions at the outer boundary!'
		!Variables derived only once - (hD_-)^8
  		vect(:, ub-3) = -vect(:,-7+ui)+8*vect(:,-6+ui) - 28*vect(:,-5+ui)+56*vect(:,-4+ui) - 70*vect(:,-3+ui)+56*vect(:,-2+ui) - 28*vect(:,-1+ui)+8*vect(:,ui)
  		vect(:, ub-2) = -vect(:,-6+ui)+8*vect(:,-5+ui) - 28*vect(:,-4+ui)+56*vect(:,-3+ui) - 70*vect(:,-2+ui)+56*vect(:,-1+ui) - 28*vect(:,ui)+8*vect(:,1+ui)
  		vect(:, ub-1) = -vect(:,-5+ui)+8*vect(:,-4+ui) - 28*vect(:,-3+ui)+56*vect(:,-2+ui) - 70*vect(:,-1+ui)+56*vect(:,ui) - 28*vect(:,1+ui)+8*vect(:,2+ui)
  		vect(:, ub-0) = -vect(:,-4+ui)+8*vect(:,-3+ui) - 28*vect(:,-2+ui)+56*vect(:,-1+ui) - 70*vect(:,ui)+56*vect(:,1+ui) - 28*vect(:,2+ui)+8*vect(:,3+ui)
  		vect(:,ub-3) = vect(:,-8+ui) - 9*vect(:,-7+ui)+36*vect(:,-6+ui) - 84*vect(:,-5+ui)+126*vect(:,-4+ui) - 126*vect(:,-3+ui)+84*vect(:,-2+ui) - 36*vect(:,-1+ui)+9*vect(:,ui)
		vect(:,ub-2) = vect(:,-7+ui) - 9*vect(:,-6+ui)+36*vect(:,-5+ui) - 84*vect(:,-4+ui)+126*vect(:,-3+ui) - 126*vect(:,-2+ui)+84*vect(:,-1+ui) - 36*vect(:,ui)+9*vect(:,1+ui)
		vect(:,ub-1) = vect(:,-6+ui) - 9*vect(:,-5+ui)+36*vect(:,-4+ui) - 84*vect(:,-3+ui)+126*vect(:,-2+ui) - 126*vect(:,-1+ui)+84*vect(:,ui) - 36*vect(:,1+ui)+9*vect(:,2+ui)
		vect(:,ub-0) = vect(:,-5+ui) - 9*vect(:,-4+ui)+36*vect(:,-3+ui) - 84*vect(:,-2+ui)+126*vect(:,-1+ui) - 126*vect(:,ui)+84*vect(:,1+ui) - 36*vect(:,2+ui)+9*vect(:,3+ui)
		!Variables derived twice - (hD_-)^9
!		vect(1:3,ub-3) = vect(1:3,-8+ui) - 9*vect(1:3,-7+ui)+36*vect(1:3,-6+ui) - 84*vect(1:3,-5+ui)+126*vect(1:3,-4+ui) - 126*vect(1:3,-3+ui)+84*vect(1:3,-2+ui) - 36*vect(1:3,-1+ui)+9*vect(1:3,ui)
!		vect(1:3,ub-2) = vect(1:3,-7+ui) - 9*vect(1:3,-6+ui)+36*vect(1:3,-5+ui) - 84*vect(1:3,-4+ui)+126*vect(1:3,-3+ui) - 126*vect(1:3,-2+ui)+84*vect(1:3,-1+ui) - 36*vect(1:3,ui)+9*vect(1:3,1+ui)
!		vect(1:3,ub-1) = vect(1:3,-6+ui) - 9*vect(1:3,-5+ui)+36*vect(1:3,-4+ui) - 84*vect(1:3,-3+ui)+126*vect(1:3,-2+ui) - 126*vect(1:3,-1+ui)+84*vect(1:3,ui) - 36*vect(1:3,1+ui)+9*vect(1:3,2+ui)
!		vect(1:3,ub-0) = vect(1:3,-5+ui) - 9*vect(1:3,-4+ui)+36*vect(1:3,-3+ui) - 84*vect(1:3,-2+ui)+126*vect(1:3,-1+ui) - 126*vect(1:3,ui)+84*vect(1:3,1+ui) - 36*vect(1:3,2+ui)+9*vect(1:3,3+ui)
!		vect(7:8,ub-3) = vect(7:8,-8+ui) - 9*vect(7:8,-7+ui)+36*vect(7:8,-6+ui) - 84*vect(7:8,-5+ui)+126*vect(7:8,-4+ui) - 126*vect(7:8,-3+ui)+84*vect(7:8,-2+ui) - 36*vect(7:8,-1+ui)+9*vect(7:8,ui)
!		vect(7:8,ub-2) = vect(7:8,-7+ui) - 9*vect(7:8,-6+ui)+36*vect(7:8,-5+ui) - 84*vect(7:8,-4+ui)+126*vect(7:8,-3+ui) - 126*vect(7:8,-2+ui)+84*vect(7:8,-1+ui) - 36*vect(7:8,ui)+9*vect(7:8,1+ui)
!		vect(7:8,ub-1) = vect(7:8,-6+ui) - 9*vect(7:8,-5+ui)+36*vect(7:8,-4+ui) - 84*vect(7:8,-3+ui)+126*vect(7:8,-2+ui) - 126*vect(7:8,-1+ui)+84*vect(7:8,ui) - 36*vect(7:8,1+ui)+9*vect(7:8,2+ui)
!		vect(7:8,ub-0) = vect(7:8,-5+ui) - 9*vect(7:8,-4+ui)+36*vect(7:8,-3+ui) - 84*vect(7:8,-2+ui)+126*vect(7:8,-1+ui) - 126*vect(7:8,ui)+84*vect(7:8,1+ui) - 36*vect(7:8,2+ui)+9*vect(7:8,3+ui)
!		vect(10,ub-3) = vect(10,-8+ui) - 9*vect(10,-7+ui)+36*vect(10,-6+ui) - 84*vect(10,-5+ui)+126*vect(10,-4+ui) - 126*vect(10,-3+ui)+84*vect(10,-2+ui) - 36*vect(10,-1+ui)+9*vect(10,ui)
!		vect(10,ub-2) = vect(10,-7+ui) - 9*vect(10,-6+ui)+36*vect(10,-5+ui) - 84*vect(10,-4+ui)+126*vect(10,-3+ui) - 126*vect(10,-2+ui)+84*vect(10,-1+ui) - 36*vect(10,ui)+9*vect(10,1+ui)
!		vect(10,ub-1) = vect(10,-6+ui) - 9*vect(10,-5+ui)+36*vect(10,-4+ui) - 84*vect(10,-3+ui)+126*vect(10,-2+ui) - 126*vect(10,-1+ui)+84*vect(10,ui) - 36*vect(10,1+ui)+9*vect(10,2+ui)
!		vect(10,ub-0) = vect(10,-5+ui) - 9*vect(10,-4+ui)+36*vect(10,-3+ui) - 84*vect(10,-2+ui)+126*vect(10,-1+ui) - 126*vect(10,ui)+84*vect(10,1+ui) - 36*vect(10,2+ui)+9*vect(10,3+ui)
	end select
	

end if

    elseif (parameters%grid%origin.eq.'stag') then
!!excision at inner boundary ! also change down 
if(parameters%grid%originbc.eq.'parity') then
	do counter = 1, abs(parameters%grid%nghost)
  		vect(:, lb + counter - 1) = + vect(:, li + abs(parameters%grid%nghost) - counter)
  		vect(5, lb + counter - 1) = - vect(5, li + abs(parameters%grid%nghost) - counter)
  		vect(7, lb + counter - 1) = - vect(7, li + abs(parameters%grid%nghost) - counter) !*!
  		vect(13:14, lb + counter - 1) = - vect(13:14, li + abs(parameters%grid%nghost) - counter)
  		vect(17, lb + counter - 1) = - vect(17, li + abs(parameters%grid%nghost) - counter)
  		!vect(n2v+1:n2v+4, lb + counter - 1) = - vect(n2v+1:n2v+4, li + abs(parameters%grid%nghost) - counter)
	enddo
endif 

if(.true.)then
! overflow boundary of 0509119
	! at xmax
!	select case (parameters%moldef%deriv_method)
	select case (parameters%grid%nghost)
! the variable's indices changed, so those derived twice need to be relabeled -- for now using the higher order outflow condition for all variables 2024/02/23
!	   case ('c2')
	   case(1)
		!Variables derived only once - (hD_-)^2
  		vect(:, ub-0) = 2*vect(:, ui-0)-vect(:, ui-1)
 		vect(:, ub-0) = 3*vect(:, ui-0)-3*vect(:, ui-1)+vect(:, ui-2) ! apparently very similar
		!Variables derived twice - (hD_-)^3
!  		vect(1:3, ub-0) = 3*vect(1:3, ui-0)-3*vect(1:3, ui-1)+vect(1:3, ui-2)
!  		vect(7:8, ub-0) = 3*vect(7:8, ui-0)-3*vect(7:8, ui-1)+vect(7:8, ui-2)
!  		vect(10, ub-0) = 3*vect(10, ui-0)-3*vect(10, ui-1)+vect(10, ui-2)

!!excision at inner boundary ! also change up and in 2nd order dissipation in numservice -- not currently updated there
if(parameters%grid%originbc.eq.'extrap')then
  		!vect(:, lb+0) = 2*vect(:, li+0)-vect(:, li+1)
  		vect(:, lb+0) = 3*vect(:, li+0)-3*vect(:, li+1)+vect(:, li+2)
endif

!	   case ('c4')
	   case(2)
		!Variables derived only once - (hD_-)^4

!  		vect(:, ub-3) = 4*vect(:, ui-2)-6*vect(:, ui-3)+4*vect(:, ui-4)-vect(:, ui-5)
 ! 		vect(:, ub-2) = 4*vect(:, ui-1)-6*vect(:, ui-2)+4*vect(:, ui-3)-vect(:, ui-4)

  		vect(:, ub-1) = 4*vect(:, ui-0)-6*vect(:, ui-1)+4*vect(:, ui-2)-vect(:, ui-3)
  		vect(:, ub-0) = 4*vect(:, ui+1)-6*vect(:, ui-0)+4*vect(:, ui-1)-vect(:, ui-2)
 		vect(:, ub-1) = 5*vect(:, ui-0)-10*vect(:, ui-1)+10*vect(:, ui-2)-5*vect(:, ui-3)+vect(:, ui-4)
 		vect(:, ub-0) = 5*vect(:, ui+1)-10*vect(:, ui-0)+10*vect(:, ui-1)-5*vect(:, ui-2)+vect(:, ui-3)
		!Variables derived twice - (hD_-)^5

  !		vect(1:3, ub-3) = 5*vect(1:3, ui-2)-10*vect(1:3, ui-3)+10*vect(1:3, ui-4)-5*vect(1:3, ui-5)+vect(1:3, ui-6)
   !    		vect(1:3, ub-2) = 5*vect(1:3, ui-1)-10*vect(1:3, ui-2)+10*vect(1:3, ui-3)-5*vect(1:3, ui-4)+vect(1:3, ui-5)
    !   		vect(7:8, ub-3) = 5*vect(7:8, ui-2)-10*vect(7:8, ui-3)+10*vect(7:8, ui-4)-5*vect(7:8, ui-5)+vect(7:8, ui-6)
!  		vect(7:8, ub-2) = 5*vect(7:8, ui-1)-10*vect(7:8, ui-2)+10*vect(7:8, ui-3)-5*vect(7:8, ui-4)+vect(7:8, ui-5)

!  		vect(1:3, ub-1) = 5*vect(1:3, ui-0)-10*vect(1:3, ui-1)+10*vect(1:3, ui-2)-5*vect(1:3, ui-3)+vect(1:3, ui-4)
!  		vect(1:3, ub-0) = 5*vect(1:3, ui+1)-10*vect(1:3, ui-0)+10*vect(1:3, ui-1)-5*vect(1:3, ui-2)+vect(1:3, ui-3)
!  		vect(7:8, ub-1) = 5*vect(7:8, ui-0)-10*vect(7:8, ui-1)+10*vect(7:8, ui-2)-5*vect(7:8, ui-3)+vect(7:8, ui-4)
!  		vect(7:8, ub-0) = 5*vect(7:8, ui+1)-10*vect(7:8, ui-0)+10*vect(7:8, ui-1)-5*vect(7:8, ui-2)+vect(7:8, ui-3)
!  		vect(10, ub-1) = 5*vect(10, ui-0)-10*vect(10, ui-1)+10*vect(10, ui-2)-5*vect(10, ui-3)+vect(10, ui-4)
!  		vect(10, ub-0) = 5*vect(10, ui+1)-10*vect(10, ui-0)+10*vect(10, ui-1)-5*vect(10, ui-2)+vect(10, ui-3)

!!excision at inner boundary ! also change up and in 4th order dissipation in numservice -- not currently updated there
if(parameters%grid%originbc.eq.'extrap')then
		!vect(:, lb+1) = 4*vect(:, li+0)-6*vect(:, li+1)+4*vect(:, li+2)-vect(:, li+3)
  		!vect(:, lb+0) = 4*vect(:, li-1)-6*vect(:, li+0)+4*vect(:, li+1)-vect(:, li+2)
  		vect(:, lb+1) = 5*vect(:, li+0)-10*vect(:, li+1)+10*vect(:, li+2)-5*vect(:, li+3)+vect(:, li+4)
  		vect(:, lb+0) = 5*vect(:, li-1)-10*vect(:, li+0)+10*vect(:, li+1)-5*vect(:, li+2)+vect(:, li+3)
endif

!	   case ('c6')
	   case(3)
!stop 'No conditions at the outer boundary!'
		!Variables derived only once - (hD_-)^6
  		vect(:, ub-2) = - vect(:,-5 + ui) + 6*vect(:,-4 + ui) - 15*vect(:,-3 + ui) + 20*vect(:,-2 + ui) - 15*vect(:,-1 + ui) + 6*vect(:,ui)
  		vect(:, ub-1) = - vect(:,-4 + ui) + 6*vect(:,-3 + ui) - 15*vect(:,-2 + ui) + 20*vect(:,-1 + ui) - 15*vect(:,ui) + 6*vect(:,1 + ui)
  		vect(:, ub-0) = - vect(:,-3 + ui) + 6*vect(:,-2 + ui) - 15*vect(:,-1 + ui) + 20*vect(:,ui) - 15*vect(:,1 + ui) + 6*vect(:,2 + ui)
  		vect(:, ub-2) = vect(:,-6 + ui) - 7*vect(:,-5 + ui) + 21*vect(:,-4 + ui) - 35*vect(:,-3 + ui) + 35*vect(:,-2 + ui) - 21*vect(:,-1 + ui) + 7*vect(:,ui)
  		vect(:, ub-1) = vect(:,-5 + ui) - 7*vect(:,-4 + ui) + 21*vect(:,-3 + ui) - 35*vect(:,-2 + ui) + 35*vect(:,-1 + ui) - 21*vect(:,ui) + 7*vect(:,1 + ui)
  		vect(:, ub-0) = vect(:,-4 + ui) - 7*vect(:,-3 + ui) + 21*vect(:,-2 + ui) - 35*vect(:,-1 + ui) + 35*vect(:,ui) - 21*vect(:,1 + ui) + 7*vect(:,2 + ui)
		!Variables derived twice - (hD_-)^7
!		vect(1:3, ub-2) = vect(1:3,-6 + ui) - 7*vect(1:3,-5 + ui) + 21*vect(1:3,-4 + ui) - 35*vect(1:3,-3 + ui) + 35*vect(1:3,-2 + ui) - 21*vect(1:3,-1 + ui) + 7*vect(1:3,ui)
!  		vect(1:3, ub-1) = vect(1:3,-5 + ui) - 7*vect(1:3,-4 + ui) + 21*vect(1:3,-3 + ui) - 35*vect(1:3,-2 + ui) + 35*vect(1:3,-1 + ui) - 21*vect(1:3,ui) + 7*vect(1:3,1 + ui)
!  		vect(1:3, ub-0) = vect(1:3,-4 + ui) - 7*vect(1:3,-3 + ui) + 21*vect(1:3,-2 + ui) - 35*vect(1:3,-1 + ui) + 35*vect(1:3,ui) - 21*vect(1:3,1 + ui) + 7*vect(1:3,2 + ui)
!		vect(7:8, ub-2) = vect(7:8,-6 + ui) - 7*vect(7:8,-5 + ui) + 21*vect(7:8,-4 + ui) - 35*vect(7:8,-3 + ui) + 35*vect(7:8,-2 + ui) - 21*vect(7:8,-1 + ui) + 7*vect(7:8,ui)
!  		vect(7:8, ub-1) = vect(7:8,-5 + ui) - 7*vect(7:8,-4 + ui) + 21*vect(7:8,-3 + ui) - 35*vect(7:8,-2 + ui) + 35*vect(7:8,-1 + ui) - 21*vect(7:8,ui) + 7*vect(7:8,1 + ui)
!  		vect(7:8, ub-0) = vect(7:8,-4 + ui) - 7*vect(7:8,-3 + ui) + 21*vect(7:8,-2 + ui) - 35*vect(7:8,-1 + ui) + 35*vect(7:8,ui) - 21*vect(7:8,1 + ui) + 7*vect(7:8,2 + ui)
!		vect(10, ub-2) = vect(10,-6 + ui) - 7*vect(10,-5 + ui) + 21*vect(10,-4 + ui) - 35*vect(10,-3 + ui) + 35*vect(10,-2 + ui) - 21*vect(10,-1 + ui) + 7*vect(10,ui)
!  		vect(10, ub-1) = vect(10,-5 + ui) - 7*vect(10,-4 + ui) + 21*vect(10,-3 + ui) - 35*vect(10,-2 + ui) + 35*vect(10,-1 + ui) - 21*vect(10,ui) + 7*vect(10,1 + ui)
!  		vect(10, ub-0) = vect(10,-4 + ui) - 7*vect(10,-3 + ui) + 21*vect(10,-2 + ui) - 35*vect(10,-1 + ui) + 35*vect(10,ui) - 21*vect(10,1 + ui) + 7*vect(10,2 + ui)

!	   case ('c8')
	   case(4)
!stop 'No conditions at the outer boundary!'
		!Variables derived only once - (hD_-)^8
  		vect(:, ub-3) = -vect(:,-7+ui)+8*vect(:,-6+ui) - 28*vect(:,-5+ui)+56*vect(:,-4+ui) - 70*vect(:,-3+ui)+56*vect(:,-2+ui) - 28*vect(:,-1+ui)+8*vect(:,ui)
  		vect(:, ub-2) = -vect(:,-6+ui)+8*vect(:,-5+ui) - 28*vect(:,-4+ui)+56*vect(:,-3+ui) - 70*vect(:,-2+ui)+56*vect(:,-1+ui) - 28*vect(:,ui)+8*vect(:,1+ui)
  		vect(:, ub-1) = -vect(:,-5+ui)+8*vect(:,-4+ui) - 28*vect(:,-3+ui)+56*vect(:,-2+ui) - 70*vect(:,-1+ui)+56*vect(:,ui) - 28*vect(:,1+ui)+8*vect(:,2+ui)
  		vect(:, ub-0) = -vect(:,-4+ui)+8*vect(:,-3+ui) - 28*vect(:,-2+ui)+56*vect(:,-1+ui) - 70*vect(:,ui)+56*vect(:,1+ui) - 28*vect(:,2+ui)+8*vect(:,3+ui)
  		vect(:,ub-3) = vect(:,-8+ui) - 9*vect(:,-7+ui)+36*vect(:,-6+ui) - 84*vect(:,-5+ui)+126*vect(:,-4+ui) - 126*vect(:,-3+ui)+84*vect(:,-2+ui) - 36*vect(:,-1+ui)+9*vect(:,ui)
		vect(:,ub-2) = vect(:,-7+ui) - 9*vect(:,-6+ui)+36*vect(:,-5+ui) - 84*vect(:,-4+ui)+126*vect(:,-3+ui) - 126*vect(:,-2+ui)+84*vect(:,-1+ui) - 36*vect(:,ui)+9*vect(:,1+ui)
		vect(:,ub-1) = vect(:,-6+ui) - 9*vect(:,-5+ui)+36*vect(:,-4+ui) - 84*vect(:,-3+ui)+126*vect(:,-2+ui) - 126*vect(:,-1+ui)+84*vect(:,ui) - 36*vect(:,1+ui)+9*vect(:,2+ui)
		vect(:,ub-0) = vect(:,-5+ui) - 9*vect(:,-4+ui)+36*vect(:,-3+ui) - 84*vect(:,-2+ui)+126*vect(:,-1+ui) - 126*vect(:,ui)+84*vect(:,1+ui) - 36*vect(:,2+ui)+9*vect(:,3+ui)
		!Variables derived twice - (hD_-)^9
!		vect(1:3,ub-3) = vect(1:3,-8+ui) - 9*vect(1:3,-7+ui)+36*vect(1:3,-6+ui) - 84*vect(1:3,-5+ui)+126*vect(1:3,-4+ui) - 126*vect(1:3,-3+ui)+84*vect(1:3,-2+ui) - 36*vect(1:3,-1+ui)+9*vect(1:3,ui)
!		vect(1:3,ub-2) = vect(1:3,-7+ui) - 9*vect(1:3,-6+ui)+36*vect(1:3,-5+ui) - 84*vect(1:3,-4+ui)+126*vect(1:3,-3+ui) - 126*vect(1:3,-2+ui)+84*vect(1:3,-1+ui) - 36*vect(1:3,ui)+9*vect(1:3,1+ui)
!		vect(1:3,ub-1) = vect(1:3,-6+ui) - 9*vect(1:3,-5+ui)+36*vect(1:3,-4+ui) - 84*vect(1:3,-3+ui)+126*vect(1:3,-2+ui) - 126*vect(1:3,-1+ui)+84*vect(1:3,ui) - 36*vect(1:3,1+ui)+9*vect(1:3,2+ui)
!		vect(1:3,ub-0) = vect(1:3,-5+ui) - 9*vect(1:3,-4+ui)+36*vect(1:3,-3+ui) - 84*vect(1:3,-2+ui)+126*vect(1:3,-1+ui) - 126*vect(1:3,ui)+84*vect(1:3,1+ui) - 36*vect(1:3,2+ui)+9*vect(1:3,3+ui)
!		vect(7:8,ub-3) = vect(7:8,-8+ui) - 9*vect(7:8,-7+ui)+36*vect(7:8,-6+ui) - 84*vect(7:8,-5+ui)+126*vect(7:8,-4+ui) - 126*vect(7:8,-3+ui)+84*vect(7:8,-2+ui) - 36*vect(7:8,-1+ui)+9*vect(7:8,ui)
!		vect(7:8,ub-2) = vect(7:8,-7+ui) - 9*vect(7:8,-6+ui)+36*vect(7:8,-5+ui) - 84*vect(7:8,-4+ui)+126*vect(7:8,-3+ui) - 126*vect(7:8,-2+ui)+84*vect(7:8,-1+ui) - 36*vect(7:8,ui)+9*vect(7:8,1+ui)
!		vect(7:8,ub-1) = vect(7:8,-6+ui) - 9*vect(7:8,-5+ui)+36*vect(7:8,-4+ui) - 84*vect(7:8,-3+ui)+126*vect(7:8,-2+ui) - 126*vect(7:8,-1+ui)+84*vect(7:8,ui) - 36*vect(7:8,1+ui)+9*vect(7:8,2+ui)
!		vect(7:8,ub-0) = vect(7:8,-5+ui) - 9*vect(7:8,-4+ui)+36*vect(7:8,-3+ui) - 84*vect(7:8,-2+ui)+126*vect(7:8,-1+ui) - 126*vect(7:8,ui)+84*vect(7:8,1+ui) - 36*vect(7:8,2+ui)+9*vect(7:8,3+ui)
!		vect(10,ub-3) = vect(10,-8+ui) - 9*vect(10,-7+ui)+36*vect(10,-6+ui) - 84*vect(10,-5+ui)+126*vect(10,-4+ui) - 126*vect(10,-3+ui)+84*vect(10,-2+ui) - 36*vect(10,-1+ui)+9*vect(10,ui)
!		vect(10,ub-2) = vect(10,-7+ui) - 9*vect(10,-6+ui)+36*vect(10,-5+ui) - 84*vect(10,-4+ui)+126*vect(10,-3+ui) - 126*vect(10,-2+ui)+84*vect(10,-1+ui) - 36*vect(10,ui)+9*vect(10,1+ui)
!		vect(10,ub-1) = vect(10,-6+ui) - 9*vect(10,-5+ui)+36*vect(10,-4+ui) - 84*vect(10,-3+ui)+126*vect(10,-2+ui) - 126*vect(10,-1+ui)+84*vect(10,ui) - 36*vect(10,1+ui)+9*vect(10,2+ui)
!		vect(10,ub-0) = vect(10,-5+ui) - 9*vect(10,-4+ui)+36*vect(10,-3+ui) - 84*vect(10,-2+ui)+126*vect(10,-1+ui) - 126*vect(10,ui)+84*vect(10,1+ui) - 36*vect(10,2+ui)+9*vect(10,3+ui)
	end select
endif

    elseif (parameters%grid%origin.eq.'misman') then !06-07-2016 misman
!!excision at inner boundary ! also change down
if(parameters%grid%originbc.eq.'parity') then
	do counter = 1, abs(parameters%grid%nghost)
  		vect(:, lb + counter - 1) = + vect(:, li + abs(parameters%grid%nghost) - counter)
  		vect(5, lb + counter - 1) = - vect(5, li + abs(parameters%grid%nghost) - counter)
  		vect(7, lb + counter - 1) = - vect(7, li + abs(parameters%grid%nghost) - counter) !*!
  		vect(13:14, lb + counter - 1) = - vect(13:14, li + abs(parameters%grid%nghost) - counter)
  		vect(17, lb + counter - 1) = - vect(17, li + abs(parameters%grid%nghost) - counter)
  		!vect(n2v+1:n2v+4, lb + counter - 1) = - vect(n2v+1:n2v+4, li + abs(parameters%grid%nghost) - counter)
	enddo
endif


if(.true.)then ! copied from the stag case - update!!!
! overflow boundary of 0509119
	! at xmax
!	select case (parameters%moldef%deriv_method)
	select case (parameters%grid%nghost)
! the variable's indices changed, so those derived twice need to be relabeled -- for now using the higher order outflow condition for all variables 2024/02/23
!	   case ('c2')
	   case(1)
		!Variables derived only once - (hD_-)^2
  		vect(:, ub-0) = 2*vect(:, ui-0)-vect(:, ui-1)
 		vect(:, ub-0) = 3*vect(:, ui-0)-3*vect(:, ui-1)+vect(:, ui-2) ! apparently very similar
		!Variables derived twice - (hD_-)^3
!  		vect(1:3, ub-0) = 3*vect(1:3, ui-0)-3*vect(1:3, ui-1)+vect(1:3, ui-2)
!  		vect(7:8, ub-0) = 3*vect(7:8, ui-0)-3*vect(7:8, ui-1)+vect(7:8, ui-2)
!  		vect(10, ub-0) = 3*vect(10, ui-0)-3*vect(10, ui-1)+vect(10, ui-2)

!!excision at inner boundary ! also change up and in 2nd order dissipation in numservice -- not currently updated there
if(parameters%grid%originbc.eq.'extrap')then
  		!vect(:, lb+0) = 2*vect(:, li+0)-vect(:, li+1)
  		vect(:, lb+0) = 3*vect(:, li+0)-3*vect(:, li+1)+vect(:, li+2)
endif

!	   case ('c4')
	   case(2)
		!Variables derived only once - (hD_-)^4

!  		vect(:, ub-3) = 4*vect(:, ui-2)-6*vect(:, ui-3)+4*vect(:, ui-4)-vect(:, ui-5)
 ! 		vect(:, ub-2) = 4*vect(:, ui-1)-6*vect(:, ui-2)+4*vect(:, ui-3)-vect(:, ui-4)

  		vect(:, ub-1) = 4*vect(:, ui-0)-6*vect(:, ui-1)+4*vect(:, ui-2)-vect(:, ui-3)
  		vect(:, ub-0) = 4*vect(:, ui+1)-6*vect(:, ui-0)+4*vect(:, ui-1)-vect(:, ui-2)
 		vect(:, ub-1) = 5*vect(:, ui-0)-10*vect(:, ui-1)+10*vect(:, ui-2)-5*vect(:, ui-3)+vect(:, ui-4)
 		vect(:, ub-0) = 5*vect(:, ui+1)-10*vect(:, ui-0)+10*vect(:, ui-1)-5*vect(:, ui-2)+vect(:, ui-3)
		!Variables derived twice - (hD_-)^5

  !		vect(1:3, ub-3) = 5*vect(1:3, ui-2)-10*vect(1:3, ui-3)+10*vect(1:3, ui-4)-5*vect(1:3, ui-5)+vect(1:3, ui-6)
   !    		vect(1:3, ub-2) = 5*vect(1:3, ui-1)-10*vect(1:3, ui-2)+10*vect(1:3, ui-3)-5*vect(1:3, ui-4)+vect(1:3, ui-5)
    !   		vect(7:8, ub-3) = 5*vect(7:8, ui-2)-10*vect(7:8, ui-3)+10*vect(7:8, ui-4)-5*vect(7:8, ui-5)+vect(7:8, ui-6)
!  		vect(7:8, ub-2) = 5*vect(7:8, ui-1)-10*vect(7:8, ui-2)+10*vect(7:8, ui-3)-5*vect(7:8, ui-4)+vect(7:8, ui-5)

!  		vect(1:3, ub-1) = 5*vect(1:3, ui-0)-10*vect(1:3, ui-1)+10*vect(1:3, ui-2)-5*vect(1:3, ui-3)+vect(1:3, ui-4)
!  		vect(1:3, ub-0) = 5*vect(1:3, ui+1)-10*vect(1:3, ui-0)+10*vect(1:3, ui-1)-5*vect(1:3, ui-2)+vect(1:3, ui-3)
!  		vect(7:8, ub-1) = 5*vect(7:8, ui-0)-10*vect(7:8, ui-1)+10*vect(7:8, ui-2)-5*vect(7:8, ui-3)+vect(7:8, ui-4)
!  		vect(7:8, ub-0) = 5*vect(7:8, ui+1)-10*vect(7:8, ui-0)+10*vect(7:8, ui-1)-5*vect(7:8, ui-2)+vect(7:8, ui-3)
!  		vect(10, ub-1) = 5*vect(10, ui-0)-10*vect(10, ui-1)+10*vect(10, ui-2)-5*vect(10, ui-3)+vect(10, ui-4)
!  		vect(10, ub-0) = 5*vect(10, ui+1)-10*vect(10, ui-0)+10*vect(10, ui-1)-5*vect(10, ui-2)+vect(10, ui-3)

!!excision at inner boundary ! also change up and in 4th order dissipation in numservice -- not currently updated there
if(parameters%grid%originbc.eq.'extrap')then
		!vect(:, lb+1) = 4*vect(:, li+0)-6*vect(:, li+1)+4*vect(:, li+2)-vect(:, li+3)
  		!vect(:, lb+0) = 4*vect(:, li-1)-6*vect(:, li+0)+4*vect(:, li+1)-vect(:, li+2)
  		vect(:, lb+1) = 5*vect(:, li+0)-10*vect(:, li+1)+10*vect(:, li+2)-5*vect(:, li+3)+vect(:, li+4)
  		vect(:, lb+0) = 5*vect(:, li-1)-10*vect(:, li+0)+10*vect(:, li+1)-5*vect(:, li+2)+vect(:, li+3)
endif

!	   case ('c6')
	   case(3)
!stop 'No conditions at the outer boundary!'
		!Variables derived only once - (hD_-)^6
  		vect(:, ub-2) = - vect(:,-5 + ui) + 6*vect(:,-4 + ui) - 15*vect(:,-3 + ui) + 20*vect(:,-2 + ui) - 15*vect(:,-1 + ui) + 6*vect(:,ui)
  		vect(:, ub-1) = - vect(:,-4 + ui) + 6*vect(:,-3 + ui) - 15*vect(:,-2 + ui) + 20*vect(:,-1 + ui) - 15*vect(:,ui) + 6*vect(:,1 + ui)
  		vect(:, ub-0) = - vect(:,-3 + ui) + 6*vect(:,-2 + ui) - 15*vect(:,-1 + ui) + 20*vect(:,ui) - 15*vect(:,1 + ui) + 6*vect(:,2 + ui)
  		vect(:, ub-2) = vect(:,-6 + ui) - 7*vect(:,-5 + ui) + 21*vect(:,-4 + ui) - 35*vect(:,-3 + ui) + 35*vect(:,-2 + ui) - 21*vect(:,-1 + ui) + 7*vect(:,ui)
  		vect(:, ub-1) = vect(:,-5 + ui) - 7*vect(:,-4 + ui) + 21*vect(:,-3 + ui) - 35*vect(:,-2 + ui) + 35*vect(:,-1 + ui) - 21*vect(:,ui) + 7*vect(:,1 + ui)
  		vect(:, ub-0) = vect(:,-4 + ui) - 7*vect(:,-3 + ui) + 21*vect(:,-2 + ui) - 35*vect(:,-1 + ui) + 35*vect(:,ui) - 21*vect(:,1 + ui) + 7*vect(:,2 + ui)
		!Variables derived twice - (hD_-)^7
!		vect(1:3, ub-2) = vect(1:3,-6 + ui) - 7*vect(1:3,-5 + ui) + 21*vect(1:3,-4 + ui) - 35*vect(1:3,-3 + ui) + 35*vect(1:3,-2 + ui) - 21*vect(1:3,-1 + ui) + 7*vect(1:3,ui)
!  		vect(1:3, ub-1) = vect(1:3,-5 + ui) - 7*vect(1:3,-4 + ui) + 21*vect(1:3,-3 + ui) - 35*vect(1:3,-2 + ui) + 35*vect(1:3,-1 + ui) - 21*vect(1:3,ui) + 7*vect(1:3,1 + ui)
!  		vect(1:3, ub-0) = vect(1:3,-4 + ui) - 7*vect(1:3,-3 + ui) + 21*vect(1:3,-2 + ui) - 35*vect(1:3,-1 + ui) + 35*vect(1:3,ui) - 21*vect(1:3,1 + ui) + 7*vect(1:3,2 + ui)
!		vect(7:8, ub-2) = vect(7:8,-6 + ui) - 7*vect(7:8,-5 + ui) + 21*vect(7:8,-4 + ui) - 35*vect(7:8,-3 + ui) + 35*vect(7:8,-2 + ui) - 21*vect(7:8,-1 + ui) + 7*vect(7:8,ui)
!  		vect(7:8, ub-1) = vect(7:8,-5 + ui) - 7*vect(7:8,-4 + ui) + 21*vect(7:8,-3 + ui) - 35*vect(7:8,-2 + ui) + 35*vect(7:8,-1 + ui) - 21*vect(7:8,ui) + 7*vect(7:8,1 + ui)
!  		vect(7:8, ub-0) = vect(7:8,-4 + ui) - 7*vect(7:8,-3 + ui) + 21*vect(7:8,-2 + ui) - 35*vect(7:8,-1 + ui) + 35*vect(7:8,ui) - 21*vect(7:8,1 + ui) + 7*vect(7:8,2 + ui)
!		vect(10, ub-2) = vect(10,-6 + ui) - 7*vect(10,-5 + ui) + 21*vect(10,-4 + ui) - 35*vect(10,-3 + ui) + 35*vect(10,-2 + ui) - 21*vect(10,-1 + ui) + 7*vect(10,ui)
!  		vect(10, ub-1) = vect(10,-5 + ui) - 7*vect(10,-4 + ui) + 21*vect(10,-3 + ui) - 35*vect(10,-2 + ui) + 35*vect(10,-1 + ui) - 21*vect(10,ui) + 7*vect(10,1 + ui)
!  		vect(10, ub-0) = vect(10,-4 + ui) - 7*vect(10,-3 + ui) + 21*vect(10,-2 + ui) - 35*vect(10,-1 + ui) + 35*vect(10,ui) - 21*vect(10,1 + ui) + 7*vect(10,2 + ui)

!	   case ('c8')
	   case(4)
!stop 'No conditions at the outer boundary!'
		!Variables derived only once - (hD_-)^8
  		vect(:, ub-3) = -vect(:,-7+ui)+8*vect(:,-6+ui) - 28*vect(:,-5+ui)+56*vect(:,-4+ui) - 70*vect(:,-3+ui)+56*vect(:,-2+ui) - 28*vect(:,-1+ui)+8*vect(:,ui)
  		vect(:, ub-2) = -vect(:,-6+ui)+8*vect(:,-5+ui) - 28*vect(:,-4+ui)+56*vect(:,-3+ui) - 70*vect(:,-2+ui)+56*vect(:,-1+ui) - 28*vect(:,ui)+8*vect(:,1+ui)
  		vect(:, ub-1) = -vect(:,-5+ui)+8*vect(:,-4+ui) - 28*vect(:,-3+ui)+56*vect(:,-2+ui) - 70*vect(:,-1+ui)+56*vect(:,ui) - 28*vect(:,1+ui)+8*vect(:,2+ui)
  		vect(:, ub-0) = -vect(:,-4+ui)+8*vect(:,-3+ui) - 28*vect(:,-2+ui)+56*vect(:,-1+ui) - 70*vect(:,ui)+56*vect(:,1+ui) - 28*vect(:,2+ui)+8*vect(:,3+ui)
  		vect(:,ub-3) = vect(:,-8+ui) - 9*vect(:,-7+ui)+36*vect(:,-6+ui) - 84*vect(:,-5+ui)+126*vect(:,-4+ui) - 126*vect(:,-3+ui)+84*vect(:,-2+ui) - 36*vect(:,-1+ui)+9*vect(:,ui)
		vect(:,ub-2) = vect(:,-7+ui) - 9*vect(:,-6+ui)+36*vect(:,-5+ui) - 84*vect(:,-4+ui)+126*vect(:,-3+ui) - 126*vect(:,-2+ui)+84*vect(:,-1+ui) - 36*vect(:,ui)+9*vect(:,1+ui)
		vect(:,ub-1) = vect(:,-6+ui) - 9*vect(:,-5+ui)+36*vect(:,-4+ui) - 84*vect(:,-3+ui)+126*vect(:,-2+ui) - 126*vect(:,-1+ui)+84*vect(:,ui) - 36*vect(:,1+ui)+9*vect(:,2+ui)
		vect(:,ub-0) = vect(:,-5+ui) - 9*vect(:,-4+ui)+36*vect(:,-3+ui) - 84*vect(:,-2+ui)+126*vect(:,-1+ui) - 126*vect(:,ui)+84*vect(:,1+ui) - 36*vect(:,2+ui)+9*vect(:,3+ui)
		!Variables derived twice - (hD_-)^9
!		vect(1:3,ub-3) = vect(1:3,-8+ui) - 9*vect(1:3,-7+ui)+36*vect(1:3,-6+ui) - 84*vect(1:3,-5+ui)+126*vect(1:3,-4+ui) - 126*vect(1:3,-3+ui)+84*vect(1:3,-2+ui) - 36*vect(1:3,-1+ui)+9*vect(1:3,ui)
!		vect(1:3,ub-2) = vect(1:3,-7+ui) - 9*vect(1:3,-6+ui)+36*vect(1:3,-5+ui) - 84*vect(1:3,-4+ui)+126*vect(1:3,-3+ui) - 126*vect(1:3,-2+ui)+84*vect(1:3,-1+ui) - 36*vect(1:3,ui)+9*vect(1:3,1+ui)
!		vect(1:3,ub-1) = vect(1:3,-6+ui) - 9*vect(1:3,-5+ui)+36*vect(1:3,-4+ui) - 84*vect(1:3,-3+ui)+126*vect(1:3,-2+ui) - 126*vect(1:3,-1+ui)+84*vect(1:3,ui) - 36*vect(1:3,1+ui)+9*vect(1:3,2+ui)
!		vect(1:3,ub-0) = vect(1:3,-5+ui) - 9*vect(1:3,-4+ui)+36*vect(1:3,-3+ui) - 84*vect(1:3,-2+ui)+126*vect(1:3,-1+ui) - 126*vect(1:3,ui)+84*vect(1:3,1+ui) - 36*vect(1:3,2+ui)+9*vect(1:3,3+ui)
!		vect(7:8,ub-3) = vect(7:8,-8+ui) - 9*vect(7:8,-7+ui)+36*vect(7:8,-6+ui) - 84*vect(7:8,-5+ui)+126*vect(7:8,-4+ui) - 126*vect(7:8,-3+ui)+84*vect(7:8,-2+ui) - 36*vect(7:8,-1+ui)+9*vect(7:8,ui)
!		vect(7:8,ub-2) = vect(7:8,-7+ui) - 9*vect(7:8,-6+ui)+36*vect(7:8,-5+ui) - 84*vect(7:8,-4+ui)+126*vect(7:8,-3+ui) - 126*vect(7:8,-2+ui)+84*vect(7:8,-1+ui) - 36*vect(7:8,ui)+9*vect(7:8,1+ui)
!		vect(7:8,ub-1) = vect(7:8,-6+ui) - 9*vect(7:8,-5+ui)+36*vect(7:8,-4+ui) - 84*vect(7:8,-3+ui)+126*vect(7:8,-2+ui) - 126*vect(7:8,-1+ui)+84*vect(7:8,ui) - 36*vect(7:8,1+ui)+9*vect(7:8,2+ui)
!		vect(7:8,ub-0) = vect(7:8,-5+ui) - 9*vect(7:8,-4+ui)+36*vect(7:8,-3+ui) - 84*vect(7:8,-2+ui)+126*vect(7:8,-1+ui) - 126*vect(7:8,ui)+84*vect(7:8,1+ui) - 36*vect(7:8,2+ui)+9*vect(7:8,3+ui)
!		vect(10,ub-3) = vect(10,-8+ui) - 9*vect(10,-7+ui)+36*vect(10,-6+ui) - 84*vect(10,-5+ui)+126*vect(10,-4+ui) - 126*vect(10,-3+ui)+84*vect(10,-2+ui) - 36*vect(10,-1+ui)+9*vect(10,ui)
!		vect(10,ub-2) = vect(10,-7+ui) - 9*vect(10,-6+ui)+36*vect(10,-5+ui) - 84*vect(10,-4+ui)+126*vect(10,-3+ui) - 126*vect(10,-2+ui)+84*vect(10,-1+ui) - 36*vect(10,ui)+9*vect(10,1+ui)
!		vect(10,ub-1) = vect(10,-6+ui) - 9*vect(10,-5+ui)+36*vect(10,-4+ui) - 84*vect(10,-3+ui)+126*vect(10,-2+ui) - 126*vect(10,-1+ui)+84*vect(10,ui) - 36*vect(10,1+ui)+9*vect(10,2+ui)
!		vect(10,ub-0) = vect(10,-5+ui) - 9*vect(10,-4+ui)+36*vect(10,-3+ui) - 84*vect(10,-2+ui)+126*vect(10,-1+ui) - 126*vect(10,ui)+84*vect(10,1+ui) - 36*vect(10,2+ui)+9*vect(10,3+ui)
	end select
endif

    else
	stop 'Is the grid staggered or not?' 
	!continue
    end if
  !		vect(4, ui-4:ub) = 0
    ! for the hyp-comp case

if(.true.) then
! the following is arranged to get the right behaviour of R (=RR, the 17th evolvar) at the origin: odd parity condition until collapse happens
if(parameters%timer%t.lt.parameters%slice%tcol) then    
      	do counter = 1, abs(parameters%grid%nghost)
  		vect(17, lb + counter - 1) = - vect(17, li + abs(parameters%grid%nghost) - counter)
!  		vect(18, lb + counter - 1) = + vect(18, li + abs(parameters%grid%nghost) - counter) ! this one does not make any visual difference
	enddo
endif
endif

    if (.false.) then
stop "omega not defined!"
	transition = exp(-(x-cuplog)**4/(x-charm)**4) ! cuplog used as ri and charm as ra (see Anil's thesis) !!!!!
	vect(8,:) = vect(8,:)*transition -x*root/(aa*aa*coefom)*(1-transition) !+ aa*vect(7,:)*coefalpha*vect(1,:)/vect(2,:)/coefgrr*domega*(1-transition)
    endif
    ! for the hyp-comp case
if(.false.)then !to avoid sobreescritura
if (parameters%moldef%stencil_boundary.ne.'right'.and.parameters%moldef%stencil_boundary.ne.'both') then
    coeff  = -1 ! dphi=coeff*ppi for the incoming wave
    if (parameters%id%id_type.eq.'gbssn-hyperb') coeff = -(sqrt(parameters%physics%aa**2+x(ui)*x(ui))-x(ui))/sqrt(parameters%physics%aa**2+x(ui)*x(ui))
  u = vect
!include '../source.f90'
	select case (parameters%moldef%deriv_method)
	   case ('c2')
	  	!vect(:, ub) = 2*vect(:,ui) - vect(:,ui-1)
 	!vect(:, ub) = vect(:,ui-1) + 2*dx*coeff*rhs(:,ui) !???????????????????????????????
		! for waveeq
 	vect(10, ub) = vect(10,ui-1) + 2*dx*coeff*vect(11,ui)
  		vect(11, ub) = 2*vect(11,ui) - vect(11,ui-1)
	   case DEFAULT
		stop 'Outer bcs not imposed for 4th, 6th and 8th order.'
	end select
end if
endif !to avoid sobreescritura

!bc
    case ('none')
  continue

!bc
    case DEFAULT
  STOP 'boundary condition unknown'

    end select
  end subroutine setboundaries



  subroutine moldeallocate
    if(tmpsteps>0) deallocate(tmp)
  end subroutine moldeallocate


end module mol
