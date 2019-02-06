! file    flrw.c
! author  Paul Lasky
! date    2014.02.20
! desc    FLRW initial data

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine FLRW_InitialData (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer   :: i, j, k
  integer, parameter :: dp = 8
  real(dp) :: a0, asq, adot, kvalue, z0
  real(dp) :: rho0, rhostar, perturb_rho0, perturb_v0
  real(dp) :: box_length, wavelength
  real(dp) :: perturb_phi, phidot, phi_offset, amp
  real(dp), parameter :: pi = 4.*atan(1.)
  real(dp) :: f, df1, df2, df3, dx, dy, dz
  real(dp) :: kx, ky, kz, modk, kval
  real(dp) :: hub, hubdot, b, L, gradH
  real(dp) :: d2chi(3,3)
  real(dp) :: cosky,cosky2,sinky,sinky2
  real(dp) :: rho0denom,rho0num

  real(dp), dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: phi_gs, delta_gs, chi_gs, rc_gs
  real(dp), dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: dxdxchi_gs,dxdychi_gs,dxdzchi_gs,dydychi_gs,dydzchi_gs,dzdzchi_gs
  real(dp), dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3),3) :: delta_vel_gs
  real(dp), dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: phi, delta, chi, rc
  real(dp), dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: dxdxchi,dxdychi,dxdzchi,dydychi,dydzchi,dzdzchi
  real(dp), dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),3) :: delta_vel

  logical   :: lapse, dtlapse, shift, data, hydro, perturb_x, perturb_y, perturb_z, perturb_all,&
       cmb_like, single_mode, perturb, synch_comov, framedrag, printy
  integer :: dr_unit, dv_unit1, dv_unit2, dv_unit3, p_unit, chi_unit, rc_unit
  integer :: xx_unit, xy_unit, xz_unit, yy_unit, yz_unit, zz_unit
  integer :: res,il,jl,kl,iu,ju,ku
  character(len=100) :: deltafile, vel1file, vel2file, vel3file, phifile, chifile, rcfile
  character(len=100) :: dxdxchifile,dxdychifile,dxdzchifile,dydychifile,dydzchifile,dzdzchifile
  character(len=100) :: dir, ics_dir, ics_type
  integer :: ics_dir_len, ics_type_len, ip1,im1,jp1,jm1,kp1,km1

  !
  ! set logicals
  !
  lapse = CCTK_EQUALS (initial_lapse, "flrw")
  dtlapse = CCTK_EQUALS (initial_dtlapse, "flrw")
  shift = CCTK_EQUALS (initial_shift, "flrw")
  data  = CCTK_EQUALS (initial_data,  "flrw")
  hydro = CCTK_EQUALS (initial_hydro, "flrw")
  perturb_x = CCTK_EQUALS (FLRW_perturb_direction, "x")
  perturb_y = CCTK_EQUALS (FLRW_perturb_direction, "y")
  perturb_z = CCTK_EQUALS (FLRW_perturb_direction, "z")
  perturb_all = CCTK_EQUALS (FLRW_perturb_direction, "all")
  single_mode = CCTK_EQUALS (FLRW_perturb_type, "single_mode")
  cmb_like = CCTK_EQUALS (FLRW_perturb_type, "CMB_like")
  synch_comov = CCTK_EQUALS (FLRW_perturb_type, "synch_comoving")
  perturb = CCTK_EQUALS (FLRW_perturb, "yes")  
  framedrag = CCTK_EQUALS (do_framedrag_test, "yes") !! must have FLRW_perturb_type = "single_mode" for this
  if (framedrag) then
       if (single_mode .eqv. .False.) then
       	  print*, 'WARNING!!! You must set FLRW_perturb_type = "single_mode" to run the framedrag test. RESET THESE and THEN RUN AGAIN!!!!'
       endif
  endif
  !
  ! set parameters
  !
  a0 = 1._dp
  asq = a0*a0
  res = int(FLRW_resolution)
  box_length = FLRW_boxlength
  dx = box_length / float(res)  !! assume xmin=0. and dx=dy=dz
  dy = dx; dz = dx
  !
  if (framedrag) then
     !
     ! set hubble, rho0 from chosen HL and FLRW_boxlength
     !
     print*, ' setting FRAMEDRAG hubble and rho_init '
     hub = framedrag_HL / box_length
     rho0 = 3._dp * hub**2 / (8._dp * pi * asq)
  else
     !
     ! set hubble from chosen FLRW_init_rho
     !
     print*, ' setting hubble and rho_init '
     rho0 = FLRW_init_rho
     hub = sqrt(8._dp * pi * rho0 * asq / 3._dp)    !! this is H (conformal) from Friedmann eqns
  endif
  hubdot = -hub**2 / 2._dp    !! can show this for conformal time (H=2/eta)
  rhostar = rho0 * a0**3      !! Conserved FLRW density
  adot = hub * a0             !! a' from H
  kvalue = -adot              !! factor outside K_ij for conformal time

  print*, 'rho0 = ',rho0

  if (perturb) then
     !
     ! set parameters only required for single mode
     !
     if (single_mode) then
     	if (framedrag) then
	    ! we just keep these names to make things simple
	    print*, ' setting FRAMEDRAG amplitudes and box_length params '
	    b = phi_perturb_amplitude
	    L = box_length
	else
	    print*, ' setting perturbation amplitude and offset params '
	    amp = phi_perturb_amplitude
     	    perturb_phi = amp * rho0	
     	    phi_offset = phi_phase_offset
        endif
     endif

     if (single_mode .or. synch_comov) then
     	print*, ' setting wavelength of perturbation and wavenumbers, modk '
     	wavelength = box_length
     	kx = 2._dp*pi/wavelength
     	ky = 2._dp*pi/wavelength
     	kz = 2._dp*pi/wavelength
     	modk = sqrt(kx**2 + ky**2 + kz**2)
     endif
     !
     ! set density, velocity amplitudes
     !
     if (single_mode) then
     	perturb_rho0 = - kx**2 / (4._dp * pi * rho0 * asq) - 2._dp
     	perturb_v0 = - sqrt(a0 / ( 6._dp * pi * rhostar ))        
     endif

     delta_vel = 0._dp
     delta = 0._dp      ! initialise perturbations to zero
     phi = 0._dp   
  endif

  if (perturb .and. cmb_like) then
     print*, ' reading in CMB_LIKE initial data ...'
     !
     ! convert CCTK_STRING "describe_ics" + "FLRW_ICs_dir" to Fortran string
     !
     call CCTK_FortranString(ics_type_len,describe_ics,ics_type)
     call CCTK_FortranString(ics_dir_len,FLRW_ICs_dir,ics_dir)
     write(dir,'(a,i3.3,a)')trim(ics_dir),res,trim(ics_type)
     !
     ! filenames to read ICs from, based on FLRW_ICs_dir + describe_ics parameters + res
     !
     write(deltafile,'(a,a)')trim(dir),'/delta.dat'
     write(vel1file,'(a,a)')trim(dir),'/vel1.dat'
     write(vel2file,'(a,a)')trim(dir),'/vel2.dat'
     write(vel3file,'(a,a)')trim(dir),'/vel3.dat'
     write(phifile,'(a,a)')trim(dir),'/phi.dat'
     !
     ! open files to read phi, delta, vel perturbs from files
     !
     open(newunit=dr_unit,file=deltafile,status='old')  ! delta rho file
     open(newunit=dv_unit1,file=vel1file,status='old')  ! delta vel file [1]
     open(newunit=dv_unit2,file=vel2file,status='old')  ! delta vel file [2]
     open(newunit=dv_unit3,file=vel3file,status='old')  ! delta vel file [3] 
     open(newunit=p_unit,file=phifile,status='old')     ! phi file
     print*,'opened IC files'
     !
     do k = 1, cctk_gsh(3)
        do j = 1, cctk_gsh(2) 
           ! loop over ROW no. (i is COLUMN no.) (i,j,k) --> (column, row, z)

           !
           ! read perturbations from CMB-like ICs files to global sized (gs) arrays
           !
           read(dr_unit,*) delta_gs(:,j,k)
           read(dv_unit1,*) delta_vel_gs(:,j,k,1)
           read(dv_unit2,*) delta_vel_gs(:,j,k,2)
           read(dv_unit3,*) delta_vel_gs(:,j,k,3)
           read(p_unit,*) phi_gs(:,j,k)
        enddo
     enddo
  endif
     
  if (perturb .and. synch_comov) then
     print*, 'reading in synch_comoving initial data...'
     !
     ! convert CCTK_STRING "describe_ics" + "FLRW_ICs_dir" to Fortran string
     !
     call CCTK_FortranString(ics_type_len,describe_ics,ics_type)
     call CCTK_FortranString(ics_dir_len,FLRW_ICs_dir,ics_dir)
     write(dir,'(a,i3.3,a)')trim(ics_dir),res,trim(ics_type)
     !print*,trim(dir)
     
     !
     ! filenames to read ICs from, based on FLRW_ICs_dir + describe_ics parameters + res
     !
     write(deltafile,'(a,a)')trim(dir),'/delta.dat'
     write(chifile,'(a,a)')trim(dir),'/chi.dat'
     write(rcfile,'(a,a)')trim(dir),'/rc.dat'

     write(dxdxchifile,'(a,a)')trim(dir),'/dxdxchi.dat'
     write(dxdychifile,'(a,a)')trim(dir),'/dxdychi.dat'
     write(dxdzchifile,'(a,a)')trim(dir),'/dxdzchi.dat'
     write(dydychifile,'(a,a)')trim(dir),'/dydychi.dat'
     write(dydzchifile,'(a,a)')trim(dir),'/dydzchi.dat'
     write(dzdzchifile,'(a,a)')trim(dir),'/dzdzchi.dat'

     !
     ! open files to read delta, chi, rc perturbs from files
     !
     open(newunit=dr_unit,file=deltafile,status='old')  ! delta rho file
     open(newunit=chi_unit,file=chifile,status='old')   ! chi file file [1]
     open(newunit=rc_unit,file=rcfile,status='old')     ! R_c file [2]
     open(newunit=xx_unit,file=dxdxchifile,status='old')! d2dx2(chi) file
     open(newunit=xy_unit,file=dxdychifile,status='old')! d2dxdy(chi) file
     open(newunit=xz_unit,file=dxdzchifile,status='old')! d2dxdx(chi) file
     open(newunit=yy_unit,file=dydychifile,status='old')! d2dy2(chi) file
     open(newunit=yz_unit,file=dydzchifile,status='old')! d2dydz(chi) file
     open(newunit=zz_unit,file=dzdzchifile,status='old')! d2dz2(chi) file

     print*,'opened IC files'

     do k = 1, cctk_gsh(3)
        do j = 1, cctk_gsh(2) 
           ! loop over ROW no. (i is COLUMN no.) (i,j,k) --> (column, row, z)

           !
           ! read perturbations from 3D files to global sized (gs) arrays
           !
           read(dr_unit,*) delta_gs(:,j,k)
           read(chi_unit,*) chi_gs(:,j,k)
           read(rc_unit,*) rc_gs(:,j,k)
	   read(xx_unit,*) dxdxchi_gs(:,j,k)
           read(xy_unit,*) dxdychi_gs(:,j,k)
           read(xz_unit,*) dxdzchi_gs(:,j,k)
           read(yy_unit,*) dydychi_gs(:,j,k)
           read(yz_unit,*) dydzchi_gs(:,j,k)
           read(zz_unit,*) dzdzchi_gs(:,j,k)
        enddo
     enddo
  endif
  
  print*, 'done reading data.'

  !
  ! indices for lower bound of local (processor) grid within global grid
  !
  il = cctk_lbnd(1) + 1
  jl = cctk_lbnd(2) + 1  ! indices output from cctk_lbnd start at 0, need to +1
  kl = cctk_lbnd(3) + 1
  !
  ! indices for upper bound of local grid within global grid
  !
  iu = cctk_ubnd(1) + 1
  ju = cctk_ubnd(2) + 1
  ku = cctk_ubnd(3) + 1
  !
  ! extract local part of global grid i.e. part this processor is using
  !
  if (perturb .and. cmb_like) then
     delta = delta_gs(il:iu, jl:ju, kl:ku)
     delta_vel = delta_vel_gs(il:iu, jl:ju, kl:ku, :)
     phi = phi_gs(il:iu, jl:ju, kl:ku)
  elseif (perturb .and. synch_comov) then
     delta = delta_gs(il:iu, jl:ju, kl:ku)
     chi = chi_gs(il:iu, jl:ju, kl:ku)
     rc = rc_gs(il:iu, jl:ju, kl:ku)
     dxdxchi = dxdxchi_gs(il:iu, jl:ju, kl:ku)
     dxdychi = dxdychi_gs(il:iu, jl:ju, kl:ku)
     dxdzchi = dxdzchi_gs(il:iu, jl:ju, kl:ku)
     dydychi = dydychi_gs(il:iu, jl:ju, kl:ku)
     dydzchi = dydzchi_gs(il:iu, jl:ju, kl:ku)
     dzdzchi = dzdzchi_gs(il:iu, jl:ju, kl:ku)
  endif


  !
  ! spatial loop
  !
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
	   printy = .False.
	   if (i==1 .and. j==1 .and. k==1) printy = .True.

           if (perturb .and. single_mode) then
              !
              ! set single mode perturbation in phi and d(phi)/dx^i
              ! use this for delta and delta_vel
              !
              if (framedrag) then
	      	 if (printy) print*, ' starting framedrag setup! '
                 ! store these cos we use them a bit in rho,vel etc
                 cosky = cos(ky * y(i,j,k))
                 cosky2 = cosky * cosky
                 sinky = sin(ky * y(i,j,k))
                 sinky2 = sinky * sinky
                 ! set grad vector \Delta_i H_j
                 ! this is the component of spatial metric \gamma_xy = \gamma_yx
                 gradH = b * cosky / (hub * L)
              else
                 if (perturb_x) then
                    f = perturb_phi * sin(kx * x(i,j,k) - phi_offset)
                    df1 = perturb_phi * kx * cos(kx * x(i,j,k) - phi_offset)
                    delta_vel(i,j,k,1) = perturb_v0 * df1
                 elseif (perturb_y) then
                    f = perturb_phi * sin(ky * y(i,j,k) - phi_offset)
                    df2 = perturb_phi * ky * cos(ky * y(i,j,k) - phi_offset)
                    delta_vel(i,j,k,2) = perturb_v0 * df2
                 elseif (perturb_z) then
                    f = perturb_phi * sin(kz * z(i,j,k) - phi_offset)
                    df3 = perturb_phi * kz * cos(kz * z(i,j,k) - phi_offset)
                    delta_vel(i,j,k,3) = perturb_v0 * df3
                 elseif (perturb_all) then
                    f = perturb_phi * (sin(kx * x(i,j,k) - phi_offset) + sin(ky * y(i,j,k) - phi_offset) + sin(kz * z(i,j,k) - phi_offset))
                    df1 = perturb_phi * kx * cos(kx * x(i,j,k) - phi_offset)
                    df2 = perturb_phi * ky * cos(ky * y(i,j,k) - phi_offset)
                    df3 = perturb_phi * kz * cos(kz * z(i,j,k) - phi_offset)
                    delta_vel(i,j,k,1) = perturb_v0 * df1
                    delta_vel(i,j,k,2) = perturb_v0 * df2
                    delta_vel(i,j,k,3) = perturb_v0 * df3
                 endif
                 phi(i,j,k) = f
                 phidot = 0._dp     !  d(phi)/dt
                 delta(i,j,k) = perturb_rho0 * f
              endif
           endif
      
           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              if (lapse) then
	      	 alp(i,j,k) = FLRW_lapse_value
                 if (perturb) then
		    !if (synch_comov .eqv. .False. .and. framedrag .eqv. .False.) then
		    if (synch_comov .eqv. .False.) then
                        if (framedrag .eqv. .False.) then
                           if (printy) print*, ' perturbing lapse! '
                           alp(i,j,k) = sqrt(1._dp + 2._dp * phi(i,j,k))
                        endif
                     endif
                  endif
               endif

              if (perturb) then
                 if (synch_comov) then
		    if (printy) print*, ' setting up synch_comov metric, ext. curvature ... '
                    ! g_ij = [1 - 2 R_c] \delta_ij + \partial_i \partial_j \chi
                    gxx(i,j,k) = asq * (1._dp - 2._dp * rc(i,j,k) + dxdxchi(i,j,k))
                    gxy(i,j,k) = asq * dxdychi(i,j,k)
                    gxz(i,j,k) = asq * dxdzchi(i,j,k)
                    gyy(i,j,k) = asq * (1._dp - 2._dp * rc(i,j,k) + dydychi(i,j,k))
                    gyz(i,j,k) = asq * dydzchi(i,j,k)
                    gzz(i,j,k) = asq * (1._dp - 2._dp * rc(i,j,k) + dzdzchi(i,j,k))

                    ! K_ij = -a' [1 - 2 R_c] \delta_ij + ( a H'/H - a' ) \partial_i \partial_j \chi
                    kxx(i,j,k) = -adot * gxx(i,j,k) + dxdxchi(i,j,k) * (a0 * hubdot / hub)
		    kxy(i,j,k) = -adot * gxy(i,j,k) + dxdychi(i,j,k) * (a0 * hubdot / hub)
		    kxz(i,j,k) = -adot * gxz(i,j,k) + dxdzchi(i,j,k) * (a0 * hubdot / hub)
		    kyy(i,j,k) = -adot * gyy(i,j,k) + dydychi(i,j,k) * (a0 * hubdot / hub)
		    kyz(i,j,k) = -adot * gyz(i,j,k) + dydzchi(i,j,k) * (a0 * hubdot / hub)
		    kzz(i,j,k) = -adot * gzz(i,j,k) + dzdzchi(i,j,k) * (a0 * hubdot / hub)
                    
                 elseif (single_mode .and. framedrag) then
		    if (printy) print*, ' setting up FRAMEDRAG metric, ext. curvature ... '
                    ! do frame-dragging potential setup - phi=psi=0, only grad vector in gamma_ij
                    gxx(i,j,k) = asq 
                    gxy(i,j,k) = gradH
                    gxz(i,j,k) = 0._dp
                    gyy(i,j,k) = asq + gradH**2
                    gyz(i,j,k) = 0._dp
                    gzz(i,j,k) = asq

                    kxx(i,j,k) = -adot
                    kxy(i,j,k) = -b * cosky / (4._dp * L)
                    kxz(i,j,k) = 0._dp
                    kyy(i,j,k) = -adot + b**2 * cosky2 / (2._dp * hub * L**2)
                    kyz(i,j,k) = 0._dp
                    kzz(i,j,k) = -adot
                 else
		    if (printy) print*, ' setting up longitudinal gauge ICs ... '
                    ! single mode or cmb-like in Longitudinal Gauge
                    gxx(i,j,k) = asq * (1._dp - 2._dp * phi(i,j,k))
                    gxy(i,j,k) = 0._dp
                    gxz(i,j,k) = 0._dp
                    gyy(i,j,k) = asq * (1._dp - 2._dp * phi(i,j,k))
                    gyz(i,j,k) = 0._dp
                    gzz(i,j,k) = asq * (1._dp - 2._dp * phi(i,j,k))
                    
                    kxx(i,j,k) = kvalue * (1._dp - 2._dp * phi(i,j,k) - phidot * a0 / adot) / alp(i,j,k)
                    kxy(i,j,k) = 0._dp
                    kxz(i,j,k) = 0._dp
                    kyy(i,j,k) = kvalue * (1._dp - 2._dp * phi(i,j,k) - phidot * a0 / adot) / alp(i,j,k)
                    kyz(i,j,k) = 0._dp
                    kzz(i,j,k) = kvalue * (1._dp - 2._dp * phi(i,j,k) - phidot * a0 / adot) / alp(i,j,k)
                 endif
              else
		if (printy) print*, ' setting up UNPERTURBED metric ... '
                 gxx(i,j,k) = asq
                 gxy(i,j,k) = 0._dp
                 gxz(i,j,k) = 0._dp
                 gyy(i,j,k) = asq
                 gyz(i,j,k) = 0._dp
                 gzz(i,j,k) = asq
                 
                 kxx(i,j,k) = kvalue
                 kxy(i,j,k) = 0._dp
                 kxz(i,j,k) = 0._dp
                 kyy(i,j,k) = kvalue
                 kyz(i,j,k) = 0._dp
                 kzz(i,j,k) = kvalue
              endif
              !
              ! time deriv of lapse -- evolution of this is specified in ADMBase.
              !
              if (dtlapse) then
                 dtalp(i,j,k) = 0._dp
              end if              

              if (shift) then
                 betax(i,j,k) = 0._dp
                 betay(i,j,k) = 0._dp
                 betaz(i,j,k) = 0._dp
              end if
              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! initialise matter to homogeneous values
                 !
		 press(i,j,k) = 0._dp ! pressure will be overwritten to P=poly_k*rho**poly_gamma in EOS_Omni
                 eps(i,j,k) = 0._dp
                 vel(i,j,k,:) = 0._dp
                 if (perturb) then
		    if (framedrag) then
		        if (printy) print*, ' setting FRAMEDRAG density, velocity ...'
                        ! set density, velocity from framedrag test (Mathematica/draft paper)
		        !
                        ! denominator and numerator in first term in notes
                        ! rho0 directly from eq.8 in draft (3.Oct.2018)
                        rho0denom = 128._dp * pi * L**2 * (16._dp * hub**2 * L**2 - 3._dp * b**2 * cosky2)
                        rho0num = (16._dp * hub**2 * L**2 - 3._dp * b**2 * cosky2)**2 - 64._dp * pi**2 * b**2 * sinky2
                        rho(i,j,k) = 3._dp * rho0num / rho0denom
			!
                        ! vel is v^i = u^i / (alp u^t) from mathematica with same u^x as draft (3.Oct.2018)
                        vel(i,j,k,1) = 8._dp * pi * b * sinky / (3._dp * b**2 * cosky2 - 16._dp * L**2 * hub**2)
			vel(i,j,k,2) = 0._dp
			vel(i,j,k,3) = 0._dp
		    else
			if (printy) print*, ' setting rho = rho_i (1+delta) ... '
			rho(i,j,k) = rho0 * ( 1._dp + delta(i,j,k) )
		    endif
		    if (framedrag .eqv. .False.) then
		       if (synch_comov .eqv. .False.) then
		           if (printy) print*, 'setting longitudinal vel=delta_vel!'
		           vel(i,j,k,:) = delta_vel(i,j,k,:)
		       endif
		    endif
                 else 
		    if (printy) print*, ' setting HOMOGENEOUS rho ... '
                    rho(i,j,k) = rho0
                 endif
		 !
              endif
              
           endif
        enddo
     enddo
  enddo
    
  
  if (CCTK_EQUALS (metric_type, "physical")) then
     ! do nothing
  else
     call CCTK_WARN (0, "Unknown value of ADMBase::metric_type -- FLRW only set-up for metric_type = physical")
  end if
  
end subroutine FLRW_InitialData



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! some functions to do some handy things !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! subroutine to return the 2nd order approx of partial 2nd deriv one variable e.g. d2/dx2(f)
! ( having this as a function didn't compile for some reason )
!
subroutine calc_deriv2(fp1,f,fm1,h,deriv2)
  integer, parameter :: dp = 8
  real(dp), intent(in) :: fp1, f, fm1  !! values of the function at i+1,i,i-1 (or j, k)
  real(dp), intent(in) :: h            !! the spacing in whatever dimension we're in
  real(dp), intent(inout) :: deriv2      !! the second deriv approx

  deriv2 = (fp1 - 2. * f + fm1) / h**2

end subroutine calc_deriv2

!
! subrouitne to return the secon *mixed* derivative, e.g. d2/dxdy(f) (general to any denom mix)
! ( compiling this as a function, calls from FLRW_InitialData didn't work either )
!
subroutine calc_deriv2_mix(f,f_xp1yp1,f_xm1ym1,f_xp1,f_xm1,f_yp1,f_ym1,dx,dy,deriv2_mix)
  integer, parameter :: dp = 8
  real(dp), intent(in) :: dx, dy ! grid spacing in the two dimensions
  real(dp), intent(in) :: f, f_xp1yp1, f_xm1ym1 ! f(x,y), f(x+dx,y+dy), f(x-dx,y-dy)
  real(dp), intent(in) :: f_xp1, f_xm1, f_yp1, f_ym1 ! f(x+dx,y), f(x-dx,y), f(x,y+dy), f(x,y-dy)
  real(dp), intent(inout) :: deriv2_mix

  real(dp) :: num, denom, d2fdx2, d2fdy2

  denom = 2. * dx * dy ! denominator

  call calc_deriv2(f_xp1,f,f_xm1,dx,d2fdx2) ! second deriv w.r.t x
  !! d2fdx2 = deriv2(f_xp1,f,f_xm1,dx) ! second deriv w.r.t x
  call calc_deriv2(f_yp1,f,f_ym1,dy,d2fdy2) ! second deriv w.r.t y
  !! d2fdy2 = deriv2(f_yp1,f,f_ym1,dy) ! second deriv w.r.t y

  num = f_xp1yp1 + f_xm1ym1 - 2. * f - dx*dx*d2fdx2 - dy*dy*d2fdy2

  deriv2_mix = num / denom  

end subroutine calc_deriv2_mix


subroutine apply_periodic(j,jp1,jm1,n)
  integer, intent(in) :: j, n ! n is size of local grid for whichever direction we're doing
  integer, intent(inout) :: jp1, jm1
  integer :: stp

  ! set (j+1), (j-1) depending on j value, implementing periodic BC's
  ! now dependent on "step" value
  if (j==1) then
     jp1 = j + 1
     jm1 = n
  elseif (j==n) then
     jp1 = 1
     jm1 = j - 1
  else
     jp1 = j + 1
     jm1 = j - 1
  endif

end subroutine apply_periodic
