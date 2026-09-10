
  null = gr ! if changed only for some rhss, also change in the dissipation implementation at the end

  chi    	= u(1,:)
  grr   	= u(2,:)
  Arr   	= u(3,:)
  K   		= u(4,:)
  Lambdar   	= u(5,:)
  alpha   	= u(6,:)
  betar   	= u(7,:)
  Theta	= u(8,:)
  rephi 	= u(9,:)
  repi 	= u(10,:)
  imphi 	= u(11,:)
  impi 	= u(12,:)
  mAr 		= u(13,:)
  mEr 		= u(14,:)
  mPhi 	= u(15,:)
  mPsi		= u(16,:)
  RR 		= u(17,:)
  TT		= u(18,:)

  du = u
  ddu = u
  u_diss = u

	last = n2v 

  if (parameters%moldef%deriv_method == 'c8') then
     call mdiff_c8(du(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
     if (abs(eps).gt.1e-16) call mdiff_diss10(u_diss(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
     call mdiff_c6(du(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
     if (abs(eps).gt.1e-16) call mdiff_diss8(u_diss(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call mdiff_c4(du(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
     if (abs(eps).gt.1e-16) call mdiff_diss6(u_diss(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
     call mdiff(du(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
     if (abs(eps).gt.1e-16) call mdiff_diss4(u_diss(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if

  dchi    	= du(1,:)
  dgrr   	= du(2,:)
  dArr   	= du(3,:)
  dK   	= du(4,:)
  dLambdar   	= du(5,:)
  dalpha   	= du(6,:)
  dbetar   	= du(7,:)
  dTheta	= du(8,:)
  drephi 	= du(9,:)
  drepi 	= du(10,:)
  dimphi 	= du(11,:)
  dimpi 	= du(12,:)
  dmAr 	= du(13,:)
  dmEr 	= du(14,:)
  dmPhi 	= du(15,:)
  dmPsi	= du(16,:)
  dRR 		= du(17,:)
  dTT		= du(18,:)

  nr0 = lbound(x,1)+parameters%grid%nghost
!  dx = x(nro+1)-x(nro)
  nro = ubound(x,1)-parameters%grid%nghost
!print*, nro, x(nro)
  dx = x(nro)-x(nro-1)


  if (parameters%moldef%spatial_order.eq.'second') then

  if (parameters%moldef%deriv_method == 'c8') then
     call m2diff_c8(ddu(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
     call m2diff_c6(ddu(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call m2diff_c4(ddu(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
     call m2diff(ddu(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if

  ddchi    	= ddu(1,:)
  ddgrr   	= ddu(2,:)
  ddArr   	= ddu(3,:)
  ddK   	= ddu(4,:)
  ddLambdar   	= ddu(5,:)
  ddalpha   	= ddu(6,:)
  ddbetar   	= ddu(7,:)
  ddTheta	= ddu(8,:)
  ddrephi 	= ddu(9,:)
  ddrepi 	= ddu(10,:)
  ddimphi 	= ddu(11,:)
  ddimpi 	= ddu(12,:)
  ddmAr 	= ddu(13,:)
  ddmEr 	= ddu(14,:)
  ddmPhi 	= ddu(15,:)
  ddmPsi	= ddu(16,:)

  	if (parameters%moldef%eqsystem.eq.'linear') then

stop 'No linear equations set.'

  	else if (parameters%moldef%eqsystem.eq.'nonlinear') then
  	
!Set unused rhs to zero
  rhs = 0; 

if(.true.)then 
  dru = u
  druin = u
  if (parameters%moldef%deriv_method == 'c8') then
  	 stop 'Off-centered derivatives for 8th order not updated.'
     call mdiffone_c8(dru(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
  	 stop 'Off-centered derivatives for 6th order not updated.'
     call mdiffone_c6(dru(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call mdiffone_c4(dru(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
     call mdiffone_c4in(druin(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
	 do counter = lbound(x,1), ubound(x,1)
	 	if (u(7,counter).gt.0) then
	 		dru(:,counter) = druin(:,counter)
	 	endif
	 enddo
  else if (parameters%moldef%deriv_method == 'c2') then
  	 stop 'Off-centered derivatives for 2nd order not updated.'
     call mdiffone(dru(1:last,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
  drchi    	= dru(1,:)
  drgrr   	= dru(2,:)
  drArr   	= dru(3,:)
  drK   	= dru(4,:)
  drLambdar   	= dru(5,:)
  dralpha   	= dru(6,:)
  drbetar   	= dru(7,:)
  drTheta	= dru(8,:)
  drrephi 	= dru(9,:)
  drrepi 	= dru(10,:)
  drimphi 	= dru(11,:)
  drimpi 	= dru(12,:)
  drmAr 	= dru(13,:)
  drmEr 	= dru(14,:)
  drmPhi 	= dru(15,:)
  drmPsi	= dru(16,:)
  drRR 	= dru(17,:)
  drTT		= dru(18,:)

!if (omt == 1) stop 'Upwind stencils for the advection terms not implemented in this system.'
!Kcmc=0._wp; omega=1._wp; domega=0._wp; ddomega=0._wp

call Vkgset(x,parameters,omega,rephi,imphi,Vkgfunc)
call Vkgprimeset(x,parameters,omega,rephi,imphi,Vkgprimefunc)


Zr = -dgrr/(2._wp*grr) + (grr*Lambdar)/2._wp + 1/x - grr**1.5_wp/x
! chi
  rhs(1,:) = null*( (-2*chi*dbetar)/3._wp + betar*dchi + (2*betar*chi*domega)/omega + (2*alpha*chi*K)/(3._wp*omega) + (2*alpha*chi*Kcmc)/(3._wp*omega) - (4*betar*chi)/(3._wp*x) + (4*alpha*chi*Theta*z4)/(3._wp*omega) )
  ! grr
  rhs(2,:) = null*( -2*alpha*Arr + betar*dgrr + (4*dbetar*grr)/3._wp - (4*betar*grr)/(3._wp*x) )
  ! Arr
  rhs(3,:) = null*( omt*betar*drArr + (1-omt)*betar*dArr + (4*Arr*dbetar)/3._wp - (2*dalpha*dchi)/3._wp - (alpha*dchi**2)/(6._wp*chi) - (2*chi*ddalpha)/3._wp + (alpha*ddchi)/3._wp + (4*alpha*chi*ddomega)/(3._wp*omega) + (Arr*betar*domega)/omega + (4*alpha*dchi*domega)/(3._wp*omega) + (7*alpha*chi*dgrr**2)/(12._wp*grr**2) - (2*alpha*Arr**2)/grr - (alpha*chi*ddgrr)/(2._wp*grr) + (chi*dalpha*dgrr)/(6._wp*grr) - (alpha*dchi*dgrr)/(12._wp*grr) - (alpha*chi*dgrr*domega)/(3._wp*omega*grr) + (2*alpha*chi*dLambdar*grr)/3._wp + (alpha*Arr*K)/omega + (alpha*Arr*Kcmc)/omega + (alpha*chi*dgrr*Lambdar)/2._wp  - (16*omega**2*alpha*chi*dimphi**2*kgf*pi)/3._wp - (16*omega**2*alpha*chi*drephi**2*kgf*pi)/3._wp - (32*omega*alpha*chi*dimphi*domega*imphi*kgf*pi)/3._wp - (16*alpha*chi*domega**2*imphi**2*kgf*pi)/3._wp + (32*omega**2*alpha*drephi*grr*imphi*kgf*mAr*pi*qp)/3._wp - (16*omega**2*alpha*grr**2*imphi**2*kgf*mAr**2*pi*qp**2)/(3._wp*chi) - (32*omega*alpha*chi*domega*drephi*kgf*pi*rephi)/3._wp - (32*omega**2*alpha*dimphi*grr*kgf*mAr*pi*qp*rephi)/3._wp - (16*alpha*chi*domega**2*kgf*pi*rephi**2)/3._wp - (16*omega**2*alpha*grr**2*kgf*mAr**2*pi*qp**2*rephi**2)/(3._wp*chi) - (2*alpha*chi)/x**2 + (2*alpha*chi*grr**1.5_wp)/x**2 - (4*Arr*betar)/(3._wp*x) + (2*chi*dalpha)/(3._wp*x) - (alpha*dchi)/(3._wp*x) - (4*alpha*chi*domega)/(3._wp*omega*x) + (2*alpha*chi*dgrr)/(3._wp*grr*x) - (5*alpha*chi*dgrr*Sqrt(grr))/(3._wp*x) - (2*alpha*chi*grr*Lambdar)/(3._wp*x) + (2*alpha*Arr*Theta*z4)/omega - (2*alpha*Arr*Theta*z4*z4c)/omega + (4*alpha*dchi*z4*Zr)/3._wp + (8*alpha*chi*domega*z4*Zr)/(3._wp*omega) )
  ! Delta K
  rhs(4,:) = null*( omt*betar*drK + (1-omt)*betar*dK + (3*omega*alpha*Arr**2)/(2._wp*grr**2) + (omega*chi*dalpha*dgrr)/grr**2 - (alpha*chi*dgrr*domega)/grr**2 + (omega*dalpha*dchi)/(2._wp*grr) - (omega*chi*ddalpha)/grr + (alpha*chi*ddomega)/grr + (3*chi*dalpha*domega)/grr - (alpha*dchi*domega)/(2._wp*grr) - (3*alpha*chi*domega**2)/(omega*grr) + (alpha*K**2)/(3._wp*omega) + (2*alpha*K*Kcmc)/(3._wp*omega) + (alpha*Kcmc**2)/(3._wp*omega)  + (8*omega**3*betar**2*dimphi**2*kgf*pi)/alpha + (8*omega**3*betar**2*drephi**2*kgf*pi)/alpha + (16*omega**2*betar**2*dimphi*domega*imphi*kgf*pi)/alpha + (8*omega*betar**2*domega**2*imphi**2*kgf*pi)/alpha - (16*omega**3*betar*dimphi*impi*kgf*pi)/alpha - (16*omega**2*betar*domega*imphi*impi*kgf*pi)/alpha + (8*omega**3*impi**2*kgf*pi)/alpha  - 16*omega**3*betar*drephi*imphi*kgf*mPhi*pi*qp + 8*omega**3*alpha*imphi**2*kgf*mPhi**2*pi*qp**2 + (16*omega**2*betar**2*domega*drephi*kgf*pi*rephi)/alpha + 16*omega**3*betar*dimphi*kgf*mPhi*pi*qp*rephi - 16*omega**3*impi*kgf*mPhi*pi*qp*rephi + (8*omega*betar**2*domega**2*kgf*pi*rephi**2)/alpha + 8*omega**3*alpha*kgf*mPhi**2*pi*qp**2*rephi**2 - (16*omega**3*betar*drephi*kgf*pi*repi)/alpha + 16*omega**3*imphi*kgf*mPhi*pi*qp*repi - (16*omega**2*betar*domega*kgf*pi*rephi*repi)/alpha + (8*omega**3*kgf*pi*repi**2)/alpha - (2*omega*chi*dalpha)/(grr*x) + (2*alpha*chi*domega)/(grr*x) + (4*alpha*K*Theta*z4)/(3._wp*omega) + (alpha*k1*Theta*z4)/omega - (alpha*k1*k2*Theta*z4)/omega + (4*alpha*Kcmc*Theta*z4)/(3._wp*omega) + (4*alpha*Theta**2*z4)/(3._wp*omega) + (2*omega*chi*dalpha*z4*z4c*Zr)/grr - (2*alpha*chi*domega*z4*z4c*Zr)/grr  -4*alpha*kgf*pi*Vkgfunc/omega )
  ! Lambdar
  rhs(5,:) = null*( omt*betar*drLambdar + (1-omt)*betar*dLambdar + (alpha*Arr*dgrr)/(2._wp*grr**3) - (betar*coeZr*dgrr**2)/(2._wp*grr**3) - (2*Arr*dalpha)/grr**2 - (3*alpha*Arr*dchi)/(chi*grr**2) - (dbetar*dgrr)/(3._wp*grr**2) - (4*alpha*Arr*domega)/(omega*grr**2) + (4*ddbetar)/(3._wp*grr) - (4*alpha*dK)/(3._wp*omega*grr) + (betar*coeZr*dgrr*Lambdar)/(2._wp*grr) - (16*omega**2*betar*dimphi**2*kgf*pi)/grr - (16*omega**2*betar*drephi**2*kgf*pi)/grr - (32*omega*betar*dimphi*domega*imphi*kgf*pi)/grr - (16*betar*domega**2*imphi**2*kgf*pi)/grr + (16*omega**2*dimphi*impi*kgf*pi)/grr + (16*omega*domega*imphi*impi*kgf*pi)/grr + (16*omega**2*betar*drephi*imphi*kgf*mAr*pi*qp)/chi + (16*omega**2*alpha*drephi*imphi*kgf*mPhi*pi*qp)/grr - (16*omega**2*alpha*imphi**2*kgf*mAr*mPhi*pi*qp**2)/chi - (32*omega*betar*domega*drephi*kgf*pi*rephi)/grr - (16*omega**2*betar*dimphi*kgf*mAr*pi*qp*rephi)/chi + (16*omega**2*impi*kgf*mAr*pi*qp*rephi)/chi - (16*omega**2*alpha*dimphi*kgf*mPhi*pi*qp*rephi)/grr - (16*betar*domega**2*kgf*pi*rephi**2)/grr - (16*omega**2*alpha*kgf*mAr*mPhi*pi*qp**2*rephi**2)/chi + (16*omega**2*drephi*kgf*pi*repi)/grr - (16*omega**2*imphi*kgf*mAr*pi*qp*repi)/chi + (16*omega*domega*kgf*pi*rephi*repi)/grr - (10*betar)/(3._wp*grr*x**2) + (2*betar*Sqrt(grr))/(3._wp*x**2) + (2*alpha*Arr)/(grr**2*x) + (4*betar*dgrr)/(3._wp*grr**2*x) + (betar*coeZr*dgrr)/(grr**2*x) + (4*dbetar)/(3._wp*grr*x) - (2*alpha*Arr)/(Sqrt(grr)*x) - (betar*coeZr*dgrr)/(Sqrt(grr)*x) + (4*dbetar*Sqrt(grr))/(3._wp*x) - (2*alpha*dTheta*z4)/(3._wp*omega*grr) + (2*alpha*domega*Theta*z4)/(omega**2*grr) - (2*dalpha*Theta*z4*z4c)/(omega*grr) - (2*alpha*k1*Zr)/(omega*grr) - (betar*coeZr*dgrr*z4*Zr)/grr**2 - (2*dbetar*z4*Zr)/(3._wp*grr) - (4*alpha*K*z4*Zr)/(3._wp*omega*grr) - (4*alpha*Kcmc*z4*Zr)/(3._wp*omega*grr) - (8*alpha*Theta*z4*Zr)/(3._wp*omega*grr) + (8*betar*z4*Zr)/(3._wp*grr*x)                          -z4*cosm*Zr/x        ) !-z4*cosm*Zr/x

!same as in 3DEin: working on it
  alphahat = sqrt((Kcmc*Kcmc*x*x/9._wp)+omega*omega)
  dalphahat = (18*omega*domega + 2*Kcmc**2*x)/(6._wp*Sqrt(9*omega**2 + Kcmc**2*x**2))
  betarhat = Kcmc*x/3._wp
  dbetarhat = Kcmc/3._wp
fwea = (2*omega*aa*xscri) !(1-x*x)**1.5_wp !(2*omega*aa*xscri)
nconst= (fwea**3*nconstv + fwea**3*nuplogv*alpha + alpha**2)  ! (fwea*(-Kcmc*xscri/3._wp)**2*nconstv + alpha**2) 
  ! alpha
!  rhs(6,:) = null*( omt*betar*dralpha + (1-omt)*betar*dalpha - (K*nconst)/omega + omega*coealphah - (alpha*betar*domega)/omega + (coealphah*Kcmc**2*x**2)/(9._wp*omega) - (omega*domega*Kcmc*x)/Sqrt(9*omega**2 + Kcmc**2*x**2) - (Kcmc**3*x**2)/(9._wp*Sqrt(9*omega**2 + Kcmc**2*x**2)) - (alpha*coealphah*Sqrt(9*omega**2 + Kcmc**2*x**2))/(3._wp*omega) + (domega*Kcmc*x*Sqrt(9*omega**2 + Kcmc**2*x**2))/(9._wp*omega) )
    rhs(6,:) = null*( -x*x*((alpha*coealphah)/omega) + x*x*(alphahat*coealphah)/omega + omt*betar*dralpha + (1-omt)*betar*dalpha - (2*alpha*betar*dalphahat)/alphahat + (alpha*betarhat*dalphahat)/alphahat + (alpha*betar*domega)/omega - (alpha*betarhat*domega)/omega - (nconst*K)/omega - (alpha**2*Kcmc)/omega - (alpha*betar**2*Kcmc)/(3._wp*omega*alphahat) + (2*alpha*betar*betarhat*Kcmc)/(3._wp*omega*alphahat) - (alpha*betarhat**2*Kcmc)/(3._wp*omega*alphahat) + (alpha**3*chi*Kcmc)/(3._wp*omega*alphahat*grr) + (2*alpha**3*chi*Sqrt(grr)*Kcmc)/(3._wp*omega*alphahat) - (2*alpha**2*Theta*z4*mu)/omega ) ! changing the coefficient of coealphah -- best convergence (for gauge waves) so far -- pbg
  ! betar
!  rhs(7,:) = null*evolshift*( omt*betar*drbetar + (1-omt)*betar*dbetar - (Kcmc**2*x)/9._wp + (3*alpha**2*chi*Lambdar)/4._wp + 3/4._wp*lambda*fwea**6*(Lambdar) + (eta*Kcmc*x)/3._wp - betar*eta ) 
    rhs(7,:) = null*evolshift*( -x*x*((betar*coealphah)/omega) + x*x*(alphahat*betar*coealphah)/omega/alpha + alphahat*dalphahat - (2*betar**2*dalphahat)/alphahat + (3*betar*betarhat*dalphahat)/alphahat - (betarhat**2*dalphahat)/alphahat + omt*betar*drbetar + (1-omt)*betar*dbetar - 2*betar*dbetarhat + betarhat*dbetarhat - (alphahat**2*domega)/omega + (betar**2*domega)/omega - (betar*betarhat*domega)/omega - (alpha*chi*dalpha)/grr + (alpha**2*dchi)/(2._wp*grr) + (3*alpha**2*chi*domega)/(omega*grr) - (2*alpha**2*chi*domega*Sqrt(grr))/omega + (2*alphahat*betar*Kcmc)/(3._wp*omega) - (betar**3*Kcmc)/(3._wp*omega*alphahat) - (2*alphahat*betarhat*Kcmc)/(3._wp*omega) + (betar**2*betarhat*Kcmc)/(omega*alphahat) - (betar*betarhat**2*Kcmc)/(omega*alphahat) + (betarhat**3*Kcmc)/(3._wp*omega*alphahat) + (alpha**2*betar*chi*Kcmc)/(3._wp*omega*alphahat*grr) - (alpha**2*betarhat*chi*Kcmc)/(3._wp*omega*alphahat*grr) + (2*alpha**2*betar*chi*Sqrt(grr)*Kcmc)/(3._wp*omega*alphahat) - (2*alpha**2*betarhat*chi*Sqrt(grr)*Kcmc)/(3._wp*omega*alphahat) + alpha**2*chi*Lambdar + fwea*lambda*Lambdar )  ! changing the coefficient of coealphah -- best convergence (for gauge waves) so far -- pbg

! Theta
  rhs(8,:) = null*z4*( omt*betar*drTheta + (1-omt)*betar*dTheta + (omega*alpha*chi*dLambdar)/2._wp + (omega*alpha*chi*dgrr**2)/(16._wp*grr**3) - (3*omega*alpha*Arr**2)/(4._wp*grr**2) - (omega*alpha*dchi*dgrr)/grr**2 - (2*alpha*chi*dgrr*domega)/grr**2 - (5*omega*alpha*dchi**2)/(4._wp*chi*grr) + (omega*alpha*ddchi)/grr + (2*alpha*chi*ddomega)/grr - (alpha*dchi*domega)/grr - (3*alpha*chi*domega**2)/(omega*grr) + (alpha*K**2)/(3._wp*omega) + (2*alpha*K*Kcmc)/(3._wp*omega) + (alpha*Kcmc**2)/(3._wp*omega)  - (4*omega**3*betar**2*dimphi**2*kgf*pi)/alpha - (4*omega**3*betar**2*drephi**2*kgf*pi)/alpha - (4*omega**3*alpha*chi*dimphi**2*kgf*pi)/grr - (4*omega**3*alpha*chi*drephi**2*kgf*pi)/grr - (8*omega**2*betar**2*dimphi*domega*imphi*kgf*pi)/alpha - (8*omega**2*alpha*chi*dimphi*domega*imphi*kgf*pi)/grr - (4*omega*betar**2*domega**2*imphi**2*kgf*pi)/alpha - (4*omega*alpha*chi*domega**2*imphi**2*kgf*pi)/grr + (8*omega**3*betar*dimphi*impi*kgf*pi)/alpha + (8*omega**2*betar*domega*imphi*impi*kgf*pi)/alpha - (4*omega**3*impi**2*kgf*pi)/alpha + 8*omega**3*alpha*drephi*imphi*kgf*mAr*pi*qp + 8*omega**3*betar*drephi*imphi*kgf*mPhi*pi*qp - (4*omega**3*alpha*grr*imphi**2*kgf*mAr**2*pi*qp**2)/chi - 4*omega**3*alpha*imphi**2*kgf*mPhi**2*pi*qp**2 - (8*omega**2*betar**2*domega*drephi*kgf*pi*rephi)/alpha - (8*omega**2*alpha*chi*domega*drephi*kgf*pi*rephi)/grr - 8*omega**3*alpha*dimphi*kgf*mAr*pi*qp*rephi - 8*omega**3*betar*dimphi*kgf*mPhi*pi*qp*rephi + 8*omega**3*impi*kgf*mPhi*pi*qp*rephi - (4*omega*betar**2*domega**2*kgf*pi*rephi**2)/alpha - (4*omega*alpha*chi*domega**2*kgf*pi*rephi**2)/grr - (4*omega**3*alpha*grr*kgf*mAr**2*pi*qp**2*rephi**2)/chi - 4*omega**3*alpha*kgf*mPhi**2*pi*qp**2*rephi**2 + (8*omega**3*betar*drephi*kgf*pi*repi)/alpha - 8*omega**3*imphi*kgf*mPhi*pi*qp*repi + (8*omega**2*betar*domega*kgf*pi*rephi*repi)/alpha - (4*omega**3*kgf*pi*repi**2)/alpha - (betar*domega*Theta)/omega + (omega*alpha*chi*dgrr)/(2._wp*grr**2*x) + (2*omega*alpha*dchi)/(grr*x) + (4*alpha*chi*domega)/(grr*x) - (omega*alpha*chi*dgrr)/(2._wp*Sqrt(grr)*x) + (omega*alpha*chi*Lambdar)/x + (4*betar*domega*Theta*z4)/omega + (4*alpha*K*Theta*z4)/(3._wp*omega) - (2*alpha*k1*Theta*z4)/omega - (alpha*k1*k2*Theta*z4)/omega + (4*alpha*Kcmc*Theta*z4)/(3._wp*omega) + (4*alpha*Theta**2*z4)/(3._wp*omega) - (3*betar*domega*Theta*z4*z4c)/omega - (alpha*K*Theta*z4*z4c)/omega - (alpha*Kcmc*Theta*z4*z4c)/omega - (2*alpha*Theta**2*z4*z4c)/omega - (omega*chi*dalpha*z4*z4c*Zr)/grr - (omega*alpha*dchi*z4*z4c*Zr)/(2._wp*grr) -4*alpha*kgf*pi*Vkgfunc/omega &
& -trafK*Theta*(1/x-1) )
! mPhi
  rhs(15,:) = 0._wp
! mPsi
  rhs(16,:) = 0._wp
 
! rephi
  rhs(9,:) = kgf*( repi )
! repi
  rhs(10,:) = kgf*( -(betar**2*ddrephi) + (betar**2*dalpha*drephi)/alpha - betar*dbetar*drephi - (3*betar**2*domega*drephi)/omega + 2*betar*drepi - (alpha**2*chi*dgrr*drephi)/grr**2 + (alpha**2*chi*ddrephi)/grr + (alpha*chi*dalpha*drephi)/grr - (alpha**2*dchi*drephi)/(2._wp*grr) - (alpha*betar*drephi*K)/omega - (alpha*betar*drephi*Kcmc)/omega - alpha**2*dmAr*imphi*qp + alpha*betar*dmPhi*imphi*qp - 2*alpha**2*dimphi*mAr*qp - alpha*dalpha*imphi*mAr*qp + (3*alpha**2*dchi*imphi*mAr*qp)/(2._wp*chi) + 2*alpha*betar*dimphi*mPhi*qp + (3*alpha*betar*domega*imphi*mPhi*qp)/omega - 2*alpha*impi*mPhi*qp + (alpha**2*imphi*K*mPhi*qp)/omega + (alpha**2*imphi*Kcmc*mPhi*qp)/omega - (betar**2*ddomega*rephi)/omega + (betar**2*dalpha*domega*rephi)/(omega*alpha) - (betar*dbetar*domega*rephi)/omega - (betar**2*domega**2*rephi)/omega**2 - (alpha**2*chi*dgrr*domega*rephi)/(omega*grr**2) + (alpha**2*chi*ddomega*rephi)/(omega*grr) + (alpha*chi*dalpha*domega*rephi)/(omega*grr) - (alpha**2*dchi*domega*rephi)/(2._wp*omega*grr) - (2*alpha**2*chi*domega**2*rephi)/(omega**2*grr) - (alpha*betar*domega*K*rephi)/omega**2 - (alpha*betar*domega*Kcmc*rephi)/omega**2  - (alpha**2*grr*mAr**2*qp**2*rephi)/chi + alpha**2*mPhi**2*qp**2*rephi - (betar*dalpha*repi)/alpha + (3*betar*domega*repi)/omega + (alpha*K*repi)/omega + (alpha*Kcmc*repi)/omega - (betar*drephi*rhs(6,:))/alpha - (betar*domega*rephi*rhs(6,:))/(omega*alpha) + (repi*rhs(6,:))/alpha + drephi*rhs(7,:) + (domega*rephi*rhs(7,:))/omega + (2*alpha**2*chi*drephi)/(grr*x) - (2*alpha**2*imphi*mAr*qp)/x + (2*alpha**2*chi*domega*rephi)/(omega*grr*x) - (2*alpha*betar*drephi*Theta*z4)/omega + (2*alpha**2*imphi*mPhi*qp*Theta*z4)/omega - (2*alpha*betar*domega*rephi*Theta*z4)/omega**2 + (2*alpha*repi*Theta*z4)/omega - alpha*imphi*qp*rhs(15,:) -((alpha**2*rephi*Vkgprimefunc)/omega**2) ) 
! imphi 
  rhs(11,:) = kgf*( impi )
! impi
  rhs(12,:) = kgf*( -(betar**2*ddimphi) + (betar**2*dalpha*dimphi)/alpha - betar*dbetar*dimphi + 2*betar*dimpi - (3*betar**2*dimphi*domega)/omega - (alpha**2*chi*dgrr*dimphi)/grr**2 + (alpha**2*chi*ddimphi)/grr + (alpha*chi*dalpha*dimphi)/grr - (alpha**2*dchi*dimphi)/(2._wp*grr) - (betar**2*ddomega*imphi)/omega + (betar**2*dalpha*domega*imphi)/(omega*alpha) - (betar*dbetar*domega*imphi)/omega - (betar**2*domega**2*imphi)/omega**2 - (alpha**2*chi*dgrr*domega*imphi)/(omega*grr**2) + (alpha**2*chi*ddomega*imphi)/(omega*grr) + (alpha*chi*dalpha*domega*imphi)/(omega*grr) - (alpha**2*dchi*domega*imphi)/(2._wp*omega*grr) - (2*alpha**2*chi*domega**2*imphi)/(omega**2*grr) - (betar*dalpha*impi)/alpha + (3*betar*domega*impi)/omega - (alpha*betar*dimphi*K)/omega - (alpha*betar*domega*imphi*K)/omega**2 + (alpha*impi*K)/omega - (alpha*betar*dimphi*Kcmc)/omega - (alpha*betar*domega*imphi*Kcmc)/omega**2 + (alpha*impi*Kcmc)/omega  + 2*alpha**2*drephi*mAr*qp - 2*alpha*betar*drephi*mPhi*qp - (alpha**2*grr*imphi*mAr**2*qp**2)/chi + alpha**2*imphi*mPhi**2*qp**2 + alpha**2*dmAr*qp*rephi - alpha*betar*dmPhi*qp*rephi + alpha*dalpha*mAr*qp*rephi - (3*alpha**2*dchi*mAr*qp*rephi)/(2._wp*chi) - (3*alpha*betar*domega*mPhi*qp*rephi)/omega - (alpha**2*K*mPhi*qp*rephi)/omega - (alpha**2*Kcmc*mPhi*qp*rephi)/omega + 2*alpha*mPhi*qp*repi - (betar*dimphi*rhs(6,:))/alpha - (betar*domega*imphi*rhs(6,:))/(omega*alpha) + (impi*rhs(6,:))/alpha + dimphi*rhs(7,:) + (domega*imphi*rhs(7,:))/omega + (2*alpha**2*chi*dimphi)/(grr*x) + (2*alpha**2*chi*domega*imphi)/(omega*grr*x) + (2*alpha**2*mAr*qp*rephi)/x - (2*alpha*betar*dimphi*Theta*z4)/omega - (2*alpha*betar*domega*imphi*Theta*z4)/omega**2 + (2*alpha*impi*Theta*z4)/omega - (2*alpha**2*mPhi*qp*rephi*Theta*z4)/omega + alpha*qp*rephi*rhs(15,:) -((alpha**2*imphi*Vkgprimefunc)/omega**2) )  

! mAr
  rhs(13,:) = 0._wp
! mEr
  rhs(14,:) = 0._wp
  
! diag variables for creating the conformal diagrams (eq A6 in 1012.3703)
  rhs(17,:) = betar*dRR + alpha*sqrt(chi/grr)*dTT 
  rhs(18,:) = betar*dTT + alpha*sqrt(chi/grr)*dRR

  constr(1,:) = (-15*chi*dgrr**2)/(8._wp*grr**3) - (3*Arr**2)/(2._wp*grr**2) + (chi*ddgrr)/grr**2 - (2*dchi*dgrr)/grr**2 - (4*chi*dgrr*domega)/(omega*grr**2) - (5*dchi**2)/(2._wp*chi*grr) + (2*ddchi)/grr + (4*chi*ddomega)/(omega*grr) - (2*dchi*domega)/(omega*grr) - (6*chi*domega**2)/(omega**2*grr) + (2*K**2)/(3._wp*omega**2) + (4*K*Kcmc)/(3._wp*omega**2) + (2*Kcmc**2)/(3._wp*omega**2)  - (8*omega**2*betar**2*dimphi**2*kgf*pi)/alpha**2 - (8*omega**2*betar**2*drephi**2*kgf*pi)/alpha**2 - (8*omega**2*chi*dimphi**2*kgf*pi)/grr - (8*omega**2*chi*drephi**2*kgf*pi)/grr - (16*omega*betar**2*dimphi*domega*imphi*kgf*pi)/alpha**2 - (16*omega*chi*dimphi*domega*imphi*kgf*pi)/grr - (8*betar**2*domega**2*imphi**2*kgf*pi)/alpha**2 - (8*chi*domega**2*imphi**2*kgf*pi)/grr + (16*omega**2*betar*dimphi*impi*kgf*pi)/alpha**2 + (16*omega*betar*domega*imphi*impi*kgf*pi)/alpha**2 - (8*omega**2*impi**2*kgf*pi)/alpha**2 + 16*omega**2*drephi*imphi*kgf*mAr*pi*qp + (16*omega**2*betar*drephi*imphi*kgf*mPhi*pi*qp)/alpha - (8*omega**2*grr*imphi**2*kgf*mAr**2*pi*qp**2)/chi - 8*omega**2*imphi**2*kgf*mPhi**2*pi*qp**2 - (16*omega*betar**2*domega*drephi*kgf*pi*rephi)/alpha**2 - (16*omega*chi*domega*drephi*kgf*pi*rephi)/grr - 16*omega**2*dimphi*kgf*mAr*pi*qp*rephi - (16*omega**2*betar*dimphi*kgf*mPhi*pi*qp*rephi)/alpha + (16*omega**2*impi*kgf*mPhi*pi*qp*rephi)/alpha - (8*betar**2*domega**2*kgf*pi*rephi**2)/alpha**2 - (8*chi*domega**2*kgf*pi*rephi**2)/grr - (8*omega**2*grr*kgf*mAr**2*pi*qp**2*rephi**2)/chi - 8*omega**2*kgf*mPhi**2*pi*qp**2*rephi**2 + (16*omega**2*betar*drephi*kgf*pi*repi)/alpha**2 - (16*omega**2*imphi*kgf*mPhi*pi*qp*repi)/alpha + (16*omega*betar*domega*kgf*pi*rephi*repi)/alpha**2 - (8*omega**2*kgf*pi*repi**2)/alpha**2 - (2*chi)/(grr*x**2) + (2*chi*Sqrt(grr))/x**2 + (5*chi*dgrr)/(grr**2*x) + (4*dchi)/(grr*x) + (8*chi*domega)/(omega*grr*x) + (8*K*Theta*z4)/(3._wp*omega**2) + (8*Kcmc*Theta*z4)/(3._wp*omega**2) + (8*Theta**2*z4)/(3._wp*omega**2)   -8*kgf*pi*Vkgfunc/omega**2 
  constr(2,:) = (-2*dK)/(3._wp*omega) - (7*Arr*dgrr)/(4._wp*grr**2) + dArr/grr - (3*Arr*dchi)/(2._wp*chi*grr) - (2*Arr*domega)/(omega*grr) - (8*omega**2*betar*dimphi**2*kgf*pi)/alpha - (8*omega**2*betar*drephi**2*kgf*pi)/alpha - (16*omega*betar*dimphi*domega*imphi*kgf*pi)/alpha - (8*betar*domega**2*imphi**2*kgf*pi)/alpha + (8*omega**2*dimphi*impi*kgf*pi)/alpha + (8*omega*domega*imphi*impi*kgf*pi)/alpha + (8*omega**2*betar*drephi*grr*imphi*kgf*mAr*pi*qp)/(alpha*chi) + 8*omega**2*drephi*imphi*kgf*mPhi*pi*qp - (8*omega**2*grr*imphi**2*kgf*mAr*mPhi*pi*qp**2)/chi - (16*omega*betar*domega*drephi*kgf*pi*rephi)/alpha - (8*omega**2*betar*dimphi*grr*kgf*mAr*pi*qp*rephi)/(alpha*chi) + (8*omega**2*grr*impi*kgf*mAr*pi*qp*rephi)/(alpha*chi) - 8*omega**2*dimphi*kgf*mPhi*pi*qp*rephi - (8*betar*domega**2*kgf*pi*rephi**2)/alpha - (8*omega**2*grr*kgf*mAr*mPhi*pi*qp**2*rephi**2)/chi + (8*omega**2*drephi*kgf*pi*repi)/alpha - (8*omega**2*grr*imphi*kgf*mAr*pi*qp*repi)/(alpha*chi) + (8*omega*domega*kgf*pi*rephi*repi)/alpha + (3*Arr)/(grr*x) - (4*dTheta*z4)/(3._wp*omega)
  constr(3,:) = -dgrr/(2._wp*grr) + (grr*Lambdar)/2._wp + 1/x - grr**1.5_wp/x

  constr(4,:) = -(alpha*dmEr) + (3*alpha*dchi*mEr)/(2._wp*chi) - 4*betar*drephi*imphi*pi*qp + 4*alpha*imphi**2*mPhi*pi*qp**2 + 4*betar*dimphi*pi*qp*rephi - 4*impi*pi*qp*rephi + 4*alpha*mPhi*pi*qp**2*rephi**2 + 4*imphi*pi*qp*repi - (2*alpha*mEr)/x

endif

! evaluation quantities newly implemented 2026/06/16: 

  ! apparhor
  constr(9,:) = grr**(-1.5_wp) - betar**2/(alpha**2*chi*Sqrt(grr)) - (dgrr*x)/(2._wp*grr**2.5_wp) - (dchi*x)/(chi*grr**1.5_wp) + (betar**2*dgrr*x)/(2._wp*alpha**2*chi*grr**1.5_wp) + (betar**2*dchi*x)/(alpha**2*chi**2*Sqrt(grr)) - (2*domega*x)/(grr**1.5_wp*omega) + (2*betar**2*domega*x)/(alpha**2*chi*Sqrt(grr)*omega) + (dgrr**2*x**2)/(16._wp*grr**3.5_wp) + (dchi*dgrr*x**2)/(4._wp*chi*grr**2.5_wp) - (betar**2*dgrr**2*x**2)/(16._wp*alpha**2*chi*grr**2.5_wp) + (dchi**2*x**2)/(4._wp*chi**2*grr**1.5_wp) - (betar**2*dchi*dgrr*x**2)/(4._wp*alpha**2*chi**2*grr**1.5_wp) - (betar**2*dchi**2*x**2)/(4._wp*alpha**2*chi**3*Sqrt(grr)) + (domega**2*x**2)/(grr**1.5_wp*omega**2) - (betar**2*domega**2*x**2)/(alpha**2*chi*Sqrt(grr)*omega**2) + (dgrr*domega*x**2)/(2._wp*grr**2.5_wp*omega) + (dchi*domega*x**2)/(chi*grr**1.5_wp*omega) - (betar**2*dgrr*domega*x**2)/(2._wp*alpha**2*chi*grr**1.5_wp*omega) - (betar**2*dchi*domega*x**2)/(alpha**2*chi**2*Sqrt(grr)*omega) - (betar*x*rhs(1,:))/(alpha**2*chi**2*Sqrt(grr)) + (betar*dgrr*x**2*rhs(1,:))/(4._wp*alpha**2*chi**2*grr**1.5_wp) + (betar*dchi*x**2*rhs(1,:))/(2._wp*alpha**2*chi**3*Sqrt(grr)) + (betar*domega*x**2*rhs(1,:))/(alpha**2*chi**2*Sqrt(grr)*omega) - (x**2*rhs(1,:)**2)/(4._wp*alpha**2*chi**3*Sqrt(grr)) - (betar*x*rhs(2,:))/(2._wp*alpha**2*chi*grr**1.5_wp) + (betar*dgrr*x**2*rhs(2,:))/(8._wp*alpha**2*chi*grr**2.5_wp) + (betar*dchi*x**2*rhs(2,:))/(4._wp*alpha**2*chi**2*grr**1.5_wp) + (betar*domega*x**2*rhs(2,:))/(2._wp*alpha**2*chi*grr**1.5_wp*omega) - (x**2*rhs(1,:)*rhs(2,:))/(4._wp*alpha**2*chi**2*grr**1.5_wp) - (x**2*rhs(2,:)**2)/(16._wp*alpha**2*chi*grr**2.5_wp)
  
  ! Misner-Sharp mass
  constr(7,:) = (x/omega/2._wp/sqrt(chi*sqrt(grr)))*( 1 - constr(9,:) )
  
  ! derivative of the Misner-Sharp mass
  if (parameters%moldef%deriv_method == 'c8') then
     call vdiff_c8(constr(7,:), constr(8,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c6') then
     call vdiff_c6(constr(7,:), constr(8,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c4') then
     call vdiff_c4(constr(7,:), constr(8,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  else if (parameters%moldef%deriv_method == 'c2') then
     call vdiff(constr(7,:), constr(8,:), x(lbound(x,1)), x(ubound(x,1)), parameters%moldef%stencil_boundary, parameters%grid%nghost)
  end if
  constr(8,:) = constr(8,:)/4._wp/pi
  
   do counter = lbound(constr,2)+parameters%grid%nghost,ubound(constr,2)-parameters%grid%nghost
   	constr(6,counter) = 0
   	do numr = lbound(constr,2)+parameters%grid%nghost, counter
   		constr(6,counter) = constr(6,counter)  + dx*kgf*4*pi*x(numr)*x(numr)*constr(5,numr)/sqrt(chi(numr)**3) !*alpha(numr)
   	end do
   end do
   
!   constr(5,:) = constr(5,:)/(x*x)
   
   constr(6,:) = constr(6,:) + mass


  ! using this slot to visualize potential function
  constr(10,:) = Vkgfunc





endif


  else if (parameters%moldef%spatial_order.eq.'first') then
!First order in space

stop "First order in space not implemented."

  else
	stop 'First or second order in space system?'

  end if


!if (mass.gt.1d-10) then
if (parameters%physics%v.gt.1d-10) then
  squ = -3*alpha*sqrt(chi/grr)/Kcmc + 0.02_wp
  rhs(1,:) = null*eps*u_diss(1,:)*squ + rhs(1,:)
  rhs(2,:) = null*eps*u_diss(2,:)*squ + rhs(2,:)
  rhs(3,:) = null*eps*u_diss(3,:)*squ + rhs(3,:)
  rhs(4,:) = null*eps*u_diss(4,:)*squ + rhs(4,:)
  rhs(5,:) = null*eps*u_diss(5,:)*squ + rhs(5,:)
  rhs(6,:) = null*eps*u_diss(6,:)*squ + rhs(6,:)
  rhs(7,:) = null*evolshift*eps*u_diss(7,:)*squ + rhs(7,:)
  rhs(8,:) = null*z4*eps*u_diss(8,:)*squ + rhs(8,:)
  rhs(9,:) = kgf*eps*u_diss(9,:)*squ + rhs(9,:)
  rhs(10,:) = kgf*eps*u_diss(10,:)*squ + rhs(10,:)
  rhs(11,:) = kgf*eps*u_diss(11,:)*squ + rhs(11,:)
  rhs(12,:) = kgf*eps*u_diss(12,:)*squ + rhs(12,:)
  rhs(13,:) = em*eps*u_diss(13,:)*squ + rhs(13,:)
  rhs(14,:) = em*eps*u_diss(14,:)*squ + rhs(14,:)
  rhs(15,:) = em*eps*u_diss(15,:)*squ + rhs(15,:)
  rhs(16,:) = em*eps*u_diss(16,:)*squ + rhs(16,:)
else
  rhs(1,:) = null*eps*u_diss(1,:) + rhs(1,:)
  rhs(2,:) = null*eps*u_diss(2,:) + rhs(2,:)
  rhs(3,:) = null*eps*u_diss(3,:) + rhs(3,:)
  rhs(4,:) = null*eps*u_diss(4,:) + rhs(4,:)
  rhs(5,:) = null*eps*u_diss(5,:) + rhs(5,:)
  rhs(6,:) = null*eps*u_diss(6,:) + rhs(6,:)
  rhs(7,:) = null*evolshift*eps*u_diss(7,:) + rhs(7,:)
  rhs(8,:) = null*z4*eps*u_diss(8,:) + rhs(8,:)
  rhs(9,:) = kgf*eps*u_diss(9,:) + rhs(9,:)
  rhs(10,:) = kgf*eps*u_diss(10,:) + rhs(10,:)
  rhs(11,:) = kgf*eps*u_diss(11,:) + rhs(11,:)
  rhs(12,:) = kgf*eps*u_diss(12,:) + rhs(12,:)
  rhs(13,:) = em*eps*u_diss(13,:) + rhs(13,:)
  rhs(14,:) = em*eps*u_diss(14,:) + rhs(14,:)
  rhs(15,:) = em*eps*u_diss(15,:) + rhs(15,:)
  rhs(16,:) = em*eps*u_diss(16,:) + rhs(16,:)
endif

