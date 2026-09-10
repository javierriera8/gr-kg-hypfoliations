program gbssnsphersym
  use mol
  use parameters_mod
  use grid, xp => gridvect
  use output
  use command_line
!  use boundaries

  implicit none

  type(parameters_type)   :: parameters, parameters_old
  integer:: i, j, k, last_step, conv_level
  character(len=120)  :: attach
  real(kind=wp) :: start, finish
  integer :: extra

  call cpu_time(start)

  call read_params(parameters, file_names, file_prefix, conv_level, attach)

if(parameters%slice%h.lt.0) then
  select case (parameters%moldef%deriv_method)
	case ('c8')
		extra=3
	case ('c6')
		extra=2
	case default
		extra=1
  end select
elseif(parameters%slice%h.gt.0) then
	extra=1
else
	stop 'Set h!'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! convergence testing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do k=1, conv_level
 parameters_old = parameters
 file_names_old = file_names
 if (conv_level == 1) then
 call change_file_names(file_names, parameters,'single')
 else
 call change_file_names(file_names, parameters, &
   & trim(attach)//achar(iachar("0")+ k))
 end if
 call main(k)
 parameters = parameters_old
 file_names = file_names_old
if(.true.) then
 parameters%grid%ncells = parameters%grid%ncells * 3.0_wp !!!!!!!!
 parameters%output%every_t_step = parameters%output%every_t_step * 3.0_wp * extra !!!!!! 9
 parameters%output%every_x_step = parameters%output%every_x_step * 3.0_wp!!!!
 parameters%timer%dt = parameters%timer%dt / 3.0_wp / extra  !!!!!!! 9
 parameters%grid%dx = parameters%grid%dx / 3.0_wp!!!!!!!!
else ! if no automatic convergence form desired
 parameters%grid%ncells = parameters%grid%ncells * 1.5_wp
 parameters%output%every_t_step = parameters%output%every_t_step * 1.5_wp * extra
 parameters%output%every_x_step = parameters%output%every_x_step
 parameters%timer%dt = parameters%timer%dt / 1.5_wp / extra
 parameters%grid%dx = parameters%grid%dx / 1.5_wp
!write(*,*) parameters%grid%ncells, parameters%timer%dt
endif
  end do

  call cpu_time(finish)

  !write(*,'(a,f8.2,a)') "Execution time = ", finish-start, " seconds"
  write(*,'(a,i2,a,i2,a,f5.2,a)') "Execution time = ", int((finish-start)/3600._wp), " h ", int(((finish-start)-int((finish-start)/3600._wp)*3600)/60._wp), " m ", dble((finish-start)-int(((finish-start)-int((finish-start)/3600._wp)*3600)/60._wp)*60)," s "

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine main(level)

implicit none

integer, intent(in) :: level
integer, parameter  :: npde = 18, nconstraints = 10 
real(kind=wp), dimension(:,:), allocatable :: u, gauge, constr
real(kind=wp), dimension(:), allocatable :: hop, psic, x
character :: delalpha*150, delbr*150, delfirst*150, delkgf*150, delz4*150, ext*4, delem*150
integer :: ns, centerindex, i, j, l, excindex
real(kind=wp) :: msmasscenter
logical :: dynupdate

real(kind=wp) ::h, k, v, lambda

ext = parameters%output%suffix

ns = 0
if (parameters%grid%origin.eq.'stag') ns = 1
if (parameters%grid%origin.eq.'misman') ns = 0 !06-07-2016 misman

allocate( u(npde, -parameters%grid%nghost : parameters%grid%ncells -ns + parameters%grid%nghost) )
allocate( constr(nconstraints, -parameters%grid%nghost : parameters%grid%ncells -ns + parameters%grid%nghost) )
allocate( gauge(6,-parameters%grid%nghost : parameters%grid%ncells -ns + parameters%grid%nghost) )
allocate( hop(-parameters%grid%nghost : parameters%grid%ncells -ns + parameters%grid%nghost) )
allocate( psic(-parameters%grid%nghost : parameters%grid%ncells -ns + parameters%grid%nghost) )

call grid_allocate(lbound(u,2), ubound(u, 2))
call new_grid(parameters)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set the gauge, initial data and the boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! include 'source.f90'

 call initialdata(u, constr, xp, gauge, parameters, psic)
! call setboundaries(u,parameters)

allocate( x(lbound(u,2) : ubound(u, 2)) )

  x = xp

!! Dynamical calculation of gauge source terms
!  centerindex = minloc(abs(x-parameters%id%center),1) ! changed to evaluate mass at the horizon
  j = lbound(constr,2) + parameters%grid%nghost
  4545 continue
  j = j+1
!  print*, constr(13,j-1), constr(13,j), constr(13,j-1)*constr(13,j)
  if (constr(nconstraints,j-1)*constr(nconstraints,j).gt.0) goto 4545
  centerindex = j
!  print*, centerindex, x(centerindex)

excindex = 0
if (abs(parameters%grid%xmin-parameters%physics%xexc).gt.1d-12) then
	excindex = minloc(abs(x-parameters%physics%xexc),1)
!	print*, "Grid index where the excision boundary will be located = ", excindex
	print*, "Excision radius = ", x(excindex)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! evolution and output loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  delalpha = 'echo Deleting output files for alpha ... & rm '//trim(file_prefix)//'*alpha*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' &'
  delbr = 'echo Deleting output files for betar and Br ... & rm '//trim(file_prefix)//'*betar*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' '//trim(file_prefix)//'*Br*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' '//trim(file_prefix)//'*derbetar*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' &'
  !delfirst = 'echo Deleting output files for the first order variables ... & rm '//trim(file_prefix)//'*der*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' &' ! change keep constraints for visualization ! undone
  delkgf = 'echo Deleting output files for phi and pi ... & rm '//trim(file_prefix)//'*phi*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' '//trim(file_prefix)//'*pi*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' '//trim(file_prefix)//'*derphi*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' &'
  delem = 'echo Deleting output files for the electric part ... & rm '//trim(file_prefix)//'*evem*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' &'
  delz4 = 'echo Deleting output files for Z4 system ... & rm '//trim(file_prefix)//'*Theta*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' &'!//' rm '//trim(file_prefix)//'*ev*Zr*_'//trim(parameters%moldef%deriv_method)//'_*'//trim(ext)//' &'
!  delalpha = trim(delalpha)
  delbr = trim(delbr)
  delfirst = trim(delfirst)
  delkgf = trim(delkgf)
  delem = trim(delem)
  delz4 = trim(delz4)

  write(*,*) 'CFL = ', parameters%timer%dt/(parameters%grid%xmax-parameters%grid%xmin)*parameters%grid%ncells
!  if (parameters%timer%dt/(parameters%grid%xmax-parameters%grid%xmin)*parameters%grid%ncells.gt.0.5d0) stop 'CFL is larger than 0.5.'

!! Dynamical calculation of gauge source terms
  msmasscenter = 0
  dynupdate = .true.
!! Dynamical calculation of gauge source terms
  last_step = ceiling((parameters%timer%end_time - parameters%timer%t) / parameters%timer%dt)
  do i = 0, last_step
 if(modulo(i, parameters%output%every_t_step) == 0) then
 call write_step(level, i, parameters, u, constr, x)!!!!!!!!!!!11
 write(*,'(F12.5)') parameters%timer%t
 end if
! call molstep(u, constr, x, gauge, parameters, psic) ! original
call molstep(u(:,excindex-parameters%grid%nghost:), constr(:,excindex-parameters%grid%nghost:), x(excindex-parameters%grid%nghost:), gauge(:,excindex-parameters%grid%nghost:), parameters, psic(excindex-parameters%grid%nghost:))
 parameters%timer%step = i+1
 if(.false.)then
!! Dynamical calculation of gauge source terms - does not really work for collapse yet, because the new source functions are implemented before the BH has formed; in the BH case the recalculation should not take place too early - when the perturbation has entered the horizon?
! if (abs(constr(11,centerindex)-msmasscenter).lt.1d-7) then !7 !12 !.and.u(7,5).gt.2d0
 if (abs(u(10,centerindex)).gt.0.3d0*parameters%id%a.and.dynupdate.eqv..true.)then
 write(*,'(a,f5.3,a,f12.8)') " Value of Misner-Sharp mass near x = ", parameters%id%center, " is: ",constr(6,minloc(abs(x-parameters%id%center),1))
 if (dynupdate) then
 parameters%physics%mass = constr(6,minloc(abs(x-parameters%id%center),1)) ! dangerous?????
 !print*, parameters%physics%mass
 call dynpsic(xp,parameters,psic)
 dynupdate = .false.
 endif
 endif
 msmasscenter = constr(6,centerindex)
!! Dynamical calculation of gauge source terms
endif


!if(.false.)then
! if (maxval(constr(1,:)) > 1E+10) then ! almost original one
! if (maxval(abs(constr(:,0:parameters%grid%ncells-ns))) > 1E+25 .or. maxval(abs(u(:,0:parameters%grid%ncells-ns))) > 1E+25) then ! most appropriate
do j=lbound(u,1), ubound(u,1)
	do l=0+1, parameters%grid%ncells-ns-1
		if(isnan(u(j,l))) then
			write(*,*) "Variable ", j, " has become NaN at x = ", x(l), " : ", u(j,l)
!if (parameters%moldef%spatial_order.eq.'second') call system(delfirst)
!if (parameters%slice%evollapse.eq.0)	call system(delalpha)
!if (parameters%slice%evolshift.eq.0)	call system(delbr)
if (parameters%physics%kgf.eq.0)	call system(delkgf)
if (parameters%physics%em.eq.0)	call system(delem)
if (parameters%physics%z4.eq.0)	call system(delz4)
  call cpu_time(finish)
  write(*,'(a,i2,a,i2,a,f5.2,a)') "Execution time = ", int((finish-start)/3600._wp), " h ", int(((finish-start)-int((finish-start)/3600._wp)*3600)/60._wp), " m ", dble((finish-start)-int(((finish-start)-int((finish-start)/3600._wp)*3600)/60._wp)*60)," s "
			stop
		endif
	enddo
enddo
 if (maxval(abs(u(:,0:parameters%grid%ncells-ns))) > 1d300) then ! normally used !1d300
! if (maxval(abs(u(1:7,0 : parameters%grid%ncells -ns))) > 5000) then ! for subsystem tests
! write(*,*) "Crash at", x(maxloc(u(2,:))) ! original
	write(*,*) "Crash at", x(maxloc(u))
!if (parameters%moldef%spatial_order.eq.'second') call system(delfirst)
!if (parameters%slice%evollapse.eq.0)	call system(delalpha)
!if (parameters%slice%evolshift.eq.0)	call system(delbr)
if (parameters%physics%kgf.eq.0)	call system(delkgf)
if (parameters%physics%em.eq.0)	call system(delem)
if (parameters%physics%z4.eq.0)	call system(delz4)

  call cpu_time(finish)
  write(*,'(a,i2,a,i2,a,f5.2,a)') "Execution time = ", int((finish-start)/3600._wp), " h ", int(((finish-start)-int((finish-start)/3600._wp)*3600)/60._wp), " m ", dble((finish-start)-int(((finish-start)-int((finish-start)/3600._wp)*3600)/60._wp)*60)," s "
 stop
 end if

!endif

  end do


!if (parameters%moldef%spatial_order.eq.'second') call system(delfirst)
!if (parameters%slice%evollapse.eq.0)	call system(delalpha)
!if (parameters%slice%evolshift.eq.0)	call system(delbr)
if (parameters%physics%kgf.eq.0)	call system(delkgf)
if (parameters%physics%em.eq.0)	call system(delem)
if (parameters%physics%z4.eq.0)	call system(delz4)

  ! Write out the last step of evolution completely to a file
!  if(parameters%id%id_type.ne.'from_data') then

 open(19, FILE=trim(file_prefix)//'last_data.txt', STATUS='REPLACE')
 write(19,*) size(u(1,:)), parameters%timer%t-parameters%timer%dt
! do i=0,parameters%grid%ncells
 do i = lbound(u,2), ubound(u,2)
 write(19,'(19(e35.24e3))') x(i), (u(j,i),j=lbound(u,1), ubound(u,1)) !! added x(i), 
 end do
 close(19)
!  end if
! '(3(e35.24e3))'
  call grid_deallocate
  deallocate(u,gauge,constr)
  call moldeallocate

end subroutine main

end program  gbssnsphersym

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The calculation of the rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcrhs(u, constr, x, gauge, rhs, parameters, psic)
  use parameters_mod
  use numservicef90
  use mol
  implicit none

  real(kind=wp), dimension(:,:), intent(inout):: u, gauge ! changed from in to inout on 24-06-2016
  real(kind=wp), dimension(:,:), intent(out):: constr
  real(kind=wp), dimension(:),intent(in):: x
  real(kind=wp), dimension(:),intent(inout):: psic
  type(parameters_type),  intent(inout) :: parameters
  real(kind=wp), dimension(size(u,1), size(u,2)), intent(out) :: rhs
!
  real(kind=wp), dimension(size(u,1), size(u,2)) :: du, ddu, u_diss, dru, druin
  real(kind=wp), dimension(size(x)) :: chi, grr, gtt, Arr, K, Lambdar, alpha, betar, Br, derchi, dergrr, dergtt, deralpha, derbetar, phikg, pikg, Theta, Zr, Lambdazr, RR, TT, rephi, repi, imphi, impi, mEr, mAr, mPhi, mPsi
  real(kind=wp), dimension(size(x)) :: dchi, dgrr, dgtt, dArr, dK, dLambdar, dalpha, dbetar, dBr, dderchi, ddergrr, ddergtt, dderalpha, dderbetar, dphikg, dpikg, dTheta, dZr, dLambdazr, dRR, dTT, drephi, drepi, dimphi, dimpi, dmEr, dmAr, dmPhi, dmPsi
  real(kind=wp), dimension(size(x)) :: drchi, drgrr, drgtt, drArr, drK, drLambdar, dralpha, drbetar, drBr, drphikg, drpikg, drTheta, drZr, drRR, drTT, drrephi, drrepi, drimphi, drimpi, drmEr, drmAr, drmPhi, drmPsi ! for advection terms
  real(kind=wp), dimension(size(x)) :: ddchi, ddgrr, ddgtt, ddArr, ddK, ddLambdar, ddalpha, ddbetar, ddBr, ddphikg, ddpikg, ddTheta, ddZr, ddrephi, ddrepi, ddimphi, ddimpi, ddmEr, ddmAr, ddmPhi, ddmPsi
  real(kind=wp), dimension(size(x)) :: auxmass, auxp, auxpp, div, squ, omega, domega, ddomega, dddomega, root, root2, coefom, ders, fwe, fwea, lcoef, Vkgfunc, Vkgprimefunc ! sqar, sqar2, sqar3, gsqarba, gsqarba2, ugsqarb, usqarb, gsqarb, ! , der -- interfering wirh dmEr -> taking this away, as I don't think it's being used
  real(kind=wp), dimension(size(x)) :: dpsic, ddpsic, coefpsi, rootp, rootp2 !, psic
  real(kind=wp), dimension(size(x)) :: dersdalpha, dersdgtt, dersdchi, dersArr, dersK, dersLambdar, dersalphau, dedersgttchi, ddersgttchi, alphahat, dalphahat, betarhat, dbetarhat
  real(kind=wp), dimension(size(x)) :: grrb, dgrrb, ddgrrb, gttb, dgttb, ddgttb, coefalpha, dcoefalpha, ddcoefalpha, coefgrr, dcoefgrr, ddcoefgrr, coefgrrb, dcoefgrrb, ddcoefgrrb, nuplog, nconst, masskgf
  real(kind=wp) :: v, eta, lambda, mu, eps, mass, cosm, cosmk, aa, cuplog, charm, xmax, k1, k2, xscri, coealphah, coealphal, coebetar, coeLambdar, coeZr, pphi, trafK, xexc, nuplogv, nconstv, k0, qp
  real(kind=wp) :: Kcmc, Ccmc, Q, rmatchmin, rmatchmax
  integer	:: evollapse, evolshift, last, nro, nr0, kgf, null, counter, numr, z4, z4c, omt, i, imatchmin, imatchmax, em, gr
  real(kind=wp) :: dx
  real(kind=wp), dimension(size(x)) :: transition, intslicing, extslicing, intshift, extshift


  integer, parameter :: n2v = 18 !change here, in the parameter file and in subroutines: initialdata, setboundaries

  cuplog = parameters%slice%cuplog
  charm = parameters%slice%charm
  evollapse = parameters%slice%evollapse
  evollapse = parameters%slice%evollapse
  evolshift = parameters%slice%evolshift
  v = parameters%physics%v
  eta = parameters%physics%eta
  lambda = parameters%physics%lambda
  mu = parameters%physics%mu
  aa = parameters%physics%aa
  mass = parameters%physics%mass
  masskgf = parameters%physics%masskgf
  cosm = parameters%physics%cosm
  gr = parameters%physics%gr
  kgf = parameters%physics%kgf
  em = parameters%physics%em
  omt = parameters%physics%omt
  cosmk = cosm/3._wp
  xmax = parameters%grid%xmax
  xscri = parameters%physics%xscri
  xexc = parameters%physics%xexc
  z4 = parameters%physics%z4
  z4c = parameters%physics%z4c
  k1 = parameters%physics%k1
  k2 = parameters%physics%k2
  coealphah = parameters%physics%coealphah
  coealphal = parameters%physics%coealphal
  nuplogv = parameters%physics%nuplogv
  nconstv = parameters%physics%nconstv
  coebetar = parameters%physics%coebetar
  coeLambdar = parameters%physics%coeLambdar
  coeZr = parameters%physics%coeZr
  pphi = parameters%physics%pphi
  trafK = parameters%physics%trafK
  Kcmc = parameters%physics%Kcmc
  Ccmc = parameters%physics%Ccmc
  Q = parameters%physics%Q
  qp = parameters%physics%qp
  rmatchmin = parameters%slice%rmatchmin
  rmatchmax = parameters%slice%rmatchmax

  if (aa.lt.0) aa = -3/Kcmc

  eps = parameters%moldef%dissipation_eps

  dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / dble(parameters%grid%ncells)
!  if (parameters%grid%origin.eq.'stag') dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / (dble(parameters%grid%ncells+1)) !change for staggered grid !RUINS CONVERGENCE!!!!!!
  if (parameters%grid%origin.eq.'misman') dx   = abs(parameters%grid%xmax - parameters%grid%xmin) / (dble(parameters%grid%ncells)+0.5_wp) !06-07-2016 misman

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

!psic=omega; Ccmc=0; mass=0

!  rootp2 = 9*Ccmc**2*psic**6+6*(Ccmc*Kcmc-3*mass)*psic**3*x**3+9*psic**2*x**4+Kcmc**2*x**6 !9*Ccmc**2 + 6*(Ccmc*Kcmc-3*mass)*x**3*psic**3 + 9*x**4*psic**4 + Kcmc**2*x**6*psic**6
do i=lbound(rootp2,1), ubound(rootp2,1) ! for some unknown reason the previous expression gives rootp2=0
  rootp2(i) = 9*Ccmc**2*psic(i)**6+6*(Ccmc*Kcmc-3*mass)*psic(i)**3*x(i)**3+9*psic(i)**2*x(i)**4+Kcmc**2*x(i)**6
enddo
do counter = lbound(u, 2), ubound(u, 2) !!crappy fix!!
	if (abs(rootp2(counter)).eq.0) rootp2(counter) = (rootp2(counter+1)+rootp2(counter-1))/2._wp
enddo
  rootp = sqrt(abs(rootp2))
!Don't use rootp or rootp2, Q=0 there.

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
  grrb = 1._wp !coefom*coefom*aa*aa/root2 !grrb
  dgrrb = 0 !aa*aa*2*coefom/root2*(-x*ddomega-coefom*(x+aa*aa*omega*domega)/root2) !dgrrb
  ddgrrb = 0 !aa*aa*2/root2*((x*ddomega+2*coefom*(x+aa*aa*omega*domega)/root2)*(x*ddomega+coefom*(x+aa*aa*omega*domega)/root2)-coefom/root2*((ddomega+x*dddomega)*root2+coefom*(1+aa*aa*domega*domega+aa*aa*omega*ddomega)-x*ddomega*(x+aa*aa*omega*domega)-2*coefom*(x+aa*aa*omega*domega)**2/root2)) !ddgrrb
!! to avoid extra noise coming from the numerical derivatives or expressions
  gttb = 1
  dgttb = 0
  ddgttb = 0
endif
!

  include './source.f90'

end subroutine calcrhs
