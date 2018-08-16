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
  real(dp) :: box_length_x, box_length_y, box_length_z, box_length, wavelength
  real(dp) :: perturb_phi, phidot, phi_offset, amp
  real(dp), parameter :: pi = 4.*atan(1.)
  real(dp) :: f, df1, df2, df3
  real(dp) :: kx, ky, kz, modk, kval
  real(dp) :: xmin, xmax
  real(dp) :: d2chi(3,3),L

  real(dp), dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: phi_gs, delta_gs, chi_gs, rc_gs
  real(dp), dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3),3) :: delta_vel_gs
  real(dp), dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: phi, delta, chi, rc
  real(dp), dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),3) :: delta_vel

  logical   :: lapse, dtlapse, shift, data, hydro, perturb_x, perturb_y, perturb_z, perturb_all,&
       cmb_like, single_mode, perturb, framedrag
  integer :: dr_unit, dv_unit1, dv_unit2, dv_unit3, p_unit, chi_unit, rc_unit
  integer :: res,il,jl,kl,iu,ju,ku
  character(len=100) :: deltafile, vel1file, vel2file, vel3file, phifile, chifile, rcfile
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
  !framedrag = CCTK_EQUALS (FLRW_perturb_type, "frame_drag_test") !! run comparison test for frame dragging potential
  
  !
  ! set parameters
  !
  a0 = 1._dp
  rho0 = FLRW_init_rho
  asq = a0*a0
  rhostar = rho0 * a0**3                     !! Conserved FLRW density
  adot = sqrt((8._dp * pi)*rho0*asq/3._dp)   !! from Friedmann eqns
  hub = adot / a0
  hubdot = -12._dp * pi * rho0 * asq / 3._dp !! H' from Friedmann eqns
  kvalue = -adot * a0                        !! check if this is really just -a' (will make no diff)
  res = int(FLRW_resolution)

  if (perturb .and. single_mode) then
     !
     ! set parameters only required for single mode
     !
     amp = phi_perturb_amplitude
     perturb_phi = amp * rho0
     phi_offset = phi_phase_offset

     wavelength = single_perturb_wavelength
     kx = 2.0_dp*pi/wavelength
     ky = 2.0_dp*pi/wavelength
     kz = 2.0_dp*pi/wavelength
     modk = sqrt(kx**2 + ky**2 + kz**2)

     if (framedrag) then
        b = amp
        L = wavelength
     endif
     !
     ! set density, velocity amplitudes
     !
     perturb_rho0 = - kx**2 / (4._dp * pi * rho0 * asq) - 2._dp
     perturb_v0 = - sqrt(a0 / ( 6._dp * pi * rhostar ))        

     delta_vel = 0._dp
     delta = 0._dp      ! initialise perturbations to zero
     phi = 0._dp   
  endif

  if (perturb .and. cmb_like) then
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

     !
     ! open files to read phi, delta, vel perturbs from files
     !
     open(newunit=dr_unit,file=deltafile,status='old')  ! delta rho file
     open(newunit=chi_unit,file=chifile,status='old')   ! chi file file [1]
     open(newunit=rc_unit,file=rcfile,status='old')     ! R_c file [2]

     print*,'opened IC files'

     do k = 1, cctk_gsh(3)
        do j = 1, cctk_gsh(2) 
           ! loop over ROW no. (i is COLUMN no.) (i,j,k) --> (column, row, z)

           !
           ! read perturbations from 3D files to global sized (gs) arrays
           !
           read(dr_unit,*) delta_gs(:,j,k)
           read(chi_unit,*) chi_gs(:,j,k,1)
           read(rc_unit,*) rc_gs(:,j,k,2)
        enddo
     enddo
  endif

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
  ! extract local part of global grid 
  !
  if (perturb .and. cmb_like) then
     delta = delta_gs(il:iu, jl:ju, kl:ku)
     delta_vel = delta_vel_gs(il:iu, jl:ju, kl:ku, :)
     phi = phi_gs(il:iu, jl:ju, kl:ku)
  elseif (perturb .and. synch_comov) then
     delta = delta_gs(il:iu, jl:ju, kl:ku)
     chi = chi_gs(il:iu, jl:ju, kl:ku)
     rc = rc_gs(il:iu, jl:ju, kl:ku)
  endif

  !
  ! spatial loop
  !
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           call apply_periodic(i,ip1,im1,nx)
           call apply_periodic(j,jp1,jm1,nx)
           call apply_periodic(k,kp1,km1,nx)

           if (perturb .and. single_mode) then
              !
              ! set single mode perturbation in phi and d(phi)/dx^i
              ! use this for delta and delta_vel
              !
              if (framedrag) then
                 ! set grad vector \Delta_i H_j, density, and vel perturb
                 gradH = (b / (hub * wavelength)) * cos(ky * y(i,j,k))
                 !delta(i,j,k) =       !!!!!!!!!!!!!!!!!! these are done in rho and vel below
                 !delta_vel(i,j,k,1) = !!!!!!!!!!!!!!!!!!
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
           ! take 2nd derivs of \chi for g_ij and K_ij
           ! ** these should be symmetric (check) **
           !
           if (synch_comov) then
              d2chi(1,1) = deriv2(chi(ip1,j,k),chi(i,j,k),chi(im1,j,k),dx)
              d2chi(1,2) = deriv2_mix(chi(i,j,k),chi(ip1,jp1,k),chi(im1,jm1,k),chi(ip1,j,k),&
                   chi(im1,j,k),chi(i,jp1,k),chi(i,jm1,k),dx,dy)
              d2chi(1,3) = deriv2_mix(chi(i,j,k),chi(ip1,j,kp1),chi(im1,j,km1),chi(ip1,j,k),&
                   chi(im1,j,k),chi(i,j,kp1),chi(i,j,km1),dx,dz)
              d2chi(2,2) = deriv2(chi(i,jp1,k),chi(i,j,k),chi(i,jm1,k),dy)
              d2chi(2,3) = deriv2_mix(chi(i,j,k),chi(i,jp1,kp1),chi(i,jm1,km1),chi(i,jp1,k),&
                   chi(i,jm1,k),chi(i,j,kp1),chi(i,j,km1),dy,dz)
              d2chi(3,3) = deriv2(chi(i,j,kp1),chi(i,j,k),chi(i,j,km1),dz)
           endif  
           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              if (lapse) then
                 if (perturb .and. synch_comov .eqv. .False. .and. framedrag .eqv. .False.) then
                    alp(i,j,k) = sqrt(1._dp + 2._dp * phi(i,j,k))
                 else
                    alp(i,j,k) = FLRW_lapse_value
                 endif
              endif

              if (perturb) then
                 if (synch_comov) then
                    ! g_ij = [1 - 2 R_c] \delta_ij + \partial_i \partial_j \chi
                    gxx(i,j,k) = asq * (1._dp - 2._dp * rc(i,j,k) + d2chi(1,1))
                    gxy(i,j,k) = asq * d2chi(1,2)
                    gxz(i,j,k) = asq * d2chi(1,3)
                    gyy(i,j,k) = asq * (1._dp - 2._dp * rc(i,j,k) + d2chi(2,2))
                    gyz(i,j,k) = asq * d2chi(2,3)
                    gzz(i,j,k) = asq * (1._dp - 2._dp * rc(i,j,k) + d2chi(3,3))

                    ! K_ij = -a' \delta_ij + H'/H \partial_i \partial_j \chi
                    kxx(i,j,k) = -adot + (hubdot * d2chi(1,1) / hub )
                    kxy(i,j,k) = (hubdot * d2chi(1,2) / hub )
                    kxz(i,j,k) = (hubdot * d2chi(1,3) / hub )
                    kyy(i,j,k) = -adot + (hubdot * d2chi(2,2) / hub )
                    kyz(i,j,k) = (hubdot * d2chi(2,3) / hub )
                    kzz(i,j,k) = -adot + (hubdot * d2chi(3,3) / hub )
                    
                 elseif (single_mode .and. framedrag) then
                    ! do frame-dragging potential setup - phi=psi=0, only grad vector in gamma_ij
                    ! g_ij = a^2 \delta_ij + \Delta_i H_j + \Delta_j H_i
                    ! second term only non-zero when \Delta_y H_x
                    gxx(i,j,k) = asq 
                    gxy(i,j,k) = asq + gradH
                    gxz(i,j,k) = 0._dp
                    gyy(i,j,k) = asq
                    gyz(i,j,k) = 0._dp
                    gzz(i,j,k) = asq

                    ! K_ij = -a' \delta_ij + H'/H \partial_i \partial_j \chi
                    kxx(i,j,k) = -adot
                    kxy(i,j,k) = -gradH * hub / (4._dp * a0) ! K_xy = -b Cos[] / 4 L a
                    kxz(i,j,k) = 0._dp
                    kyy(i,j,k) = -adot
                    kyz(i,j,k) = 0._dp
                    kzz(i,j,k) = -adot
                 else
                    ! single mode or cmb-like in Longitudinal Gauge
                    gxx(i,j,k) = asq * (1._dp - 2._dp * phi(i,j,k))
                    gxy(i,j,k) = 0._dp
                    gxz(i,j,k) = 0._dp
                    gyy(i,j,k) = asq * (1._dp - 2._dp * phi(i,j,k))
                    gyz(i,j,k) = 0._dp
                    gzz(i,j,k) = asq * (1._dp - 2._dp * phi(i,j,k))

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!!!!!!!! check K_ij b/g is -a*a' in conformal gauge... might just be -a'
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
                    kxx(i,j,k) = kvalue * (1._dp - 2._dp * phi(i,j,k) - phidot * a0 / adot) / alp(i,j,k)
                    kxy(i,j,k) = 0._dp
                    kxz(i,j,k) = 0._dp
                    kyy(i,j,k) = kvalue * (1._dp - 2._dp * phi(i,j,k) - phidot * a0 / adot) / alp(i,j,k)
                    kyz(i,j,k) = 0._dp
                    kzz(i,j,k) = kvalue * (1._dp - 2._dp * phi(i,j,k) - phidot * a0 / adot) / alp(i,j,k)
                 endif
              else
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
                 if (perturb .and. framedrag) then
                    ! set density, velocity from framedrag test (Mathematica)
                    L = wavelength
                    ! this rho0 (rest-frame) from Mathematica nb // FortranForm + vars re-named
                    rho(i,j,k) = (-18. * b**2 * L**4 * pi * a0**9 * hub**6 * sin(ky*y(i,j,k))**2)/&
                         ((-(b**2 * cos(ky*y(i,j,k))**2) + L**2 * a0**4 * hub**2)**3 * &
                         ((-32. * b**2 * cos(ky*y*i,j,k))**2)/a0 + 12. * L**2 * a0**3 * hub**4 * adot**2 - &
                         b**2 * a0 * cos(ky*y(i,j,k))**2 * hdot**2)) + &
                         ((-32. * b**2 * cos(ky*y(i,j,k))**2)/a0 + 12. * L**2 * a0**3 * hub**4 * adot**2 - &
                         b**2 * a0 * cos(ky*y(i,j,k))**2 * hdot**2) / &
                         (32. * pi * a0**3 * hub**2 * (-(b**2 * cos(ky*y(i,j,k))**2) + L**2 * a0**4 * hub**2))
                    !
                    vel(i,j,k,1) = 2._dp * b * pi * sin(ky*y(i,j,k)) / (L**2 * adot**3) !! this is up to b^3 accurate - exact form is super ugly...
                 else
                    rho(i,j,k) = rho0
                 endif
                 press(i,j,k) = 0._dp ! pressure will be overwritten to P=poly_k*rho**poly_gamma in EOS_Omni
                 eps(i,j,k) = 0._dp
                 vel(i,j,k,:) = 0._dp
                 if (perturb .and. framedrag .eqv. .False.) then
                    rho(i,j,k) = rho(i,j,k) * ( 1._dp + delta(i,j,k) )
                    if (synch_comov .eqv. .False.) then
                       vel(i,j,k,:) = delta_vel(i,j,k,:)
                    endif
                 endif
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
! function to return the 2nd order approx of partial 2nd deriv one variable e.g. d2/dx2(f)
! 
real(c_double) function deriv2(fp1,f,fm1,h)
  real(c_double) :: fp1, f, fm1  !! values of the function at i+1,i,i-1 (or j, k)
  real(c_double) :: h !! the spacing in whatever dimension we're in

  deriv2 = (fp1 - 2. * f + fm1) / h**2

end function deriv2

!
! a function to return the secon *mixed* derivative, e.g. d2/dxdy(f) (general to any denom mix)
!
real(c_double) function deriv2_mix(f,f_xp1yp1,f_xm1ym1,f_xp1,f_xm1,f_yp1,f_ym1,dx,dy)
  real(c_double), intent(in) :: dx, dy ! grid spacing in the two dimensions
  real(c_double), intent(in) :: f, f_xp1yp1, f_xm1ym1 ! f(x,y), f(x+dx,y+dy), f(x-dx,y-dy)
  real(c_double), intent(in) :: f_xp1, f_xm1, f_yp1, f_ym1 ! f(x+dx,y), f(x-dx,y), f(x,y+dy), f(x,y-dy)

  real(c_double) :: num, denom, d2fdx2, d2fdy2

  denom = 2. * dx * dy ! denominator

  d2fdx2 = deriv2(f_xp1,f,f_xm1,dx) ! second deriv w.r.t x
  d2fdy2 = deriv2(f_yp1,f,f_ym1,dy) ! second deriv w.r.t y

  num = f_xp1yp1 + f_xm1ym1 - 2. * f - dx*dx*d2fdx2 - dy*dy*d2fdy2

  deriv2_mix = num / denom  
end function deriv2_mix


subroutine apply_periodic(j,jp1,jm1,nx)
  integer, intent(in) :: j, nx
  integer, intent(inout) :: jp1, jm1
  integer :: stp

  ! set (j+1), (j-1) depending on j value, implementing periodic BC's
  ! now dependent on "step" value
  if (j==1) then
     jp1 = j + 1
     jm1 = nx
  elseif (j==nx) then
     jp1 = 1
     jm1 = j - 1
  else
     jp1 = j + 1
     jm1 = j - 1
  endif

end subroutine apply_periodic
