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
  real(dp), parameter :: pi = 3.14159265358979323846264338327
  real(dp) :: f, df1, df2, df3
  real(dp) :: kx, ky, kz, modk, kval

  real(dp), dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: phi, delta
  real(dp), dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),3) :: delta_vel
  real(dp), dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: phi_gs, delta_gs
  real(dp), dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3),3) :: delta_vel_gs

  logical   :: lapse, dtlapse, shift, data, hydro, perturb_x, perturb_y, perturb_z, perturb_all, cmb_like, single_mode
  integer :: dr_unit, dv_unit1, dv_unit2, dv_unit3, p_unit
  integer :: res, il, jl, kl, iu, ju, ku
  character(len=40) :: deltafile, vel1file, vel2file, vel3file, phifile

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

  !
  ! set parameters
  !
  !z0 = FLRW_initial_redshift
  !a0 =1._dp / (1._dp + z0)     !! want to include redshift in future. for now cannot simulate other then a0=1
  a0 = 1._dp
  rho0 = FLRW_init_rho
  asq = a0*a0
  rhostar = rho0 * asq*a0   !! Conserved FLRW density
  adot = sqrt((8._dp * pi)*rho0*asq/3._dp) !! from Friedmann eqns
  kvalue = -adot * a0
  res = int(FLRW_resolution)

  if (single_mode .and. FLRW_perturb) then

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

     ! set density, velocity amplitudes
     if (perturb_x) then
        kval = kx
     elseif (perturb_y) then
        kval = ky
     elseif (perturb_z) then
        kval = kz
     elseif (perturb_all) then
        kval = modk
     endif
!     perturb_rho0 = kval**2 / (4._dp * pi * rho0 * asq) - 2._dp
     perturb_rho0 = - kx**2 / (4._dp * pi * rho0 * asq) - 2._dp  !! testing - why kx only here?
     perturb_v0 = - sqrt(a0 / ( 6._dp * pi * rhostar ))        

     delta_vel = 0._dp
     delta = 0._dp      ! initialise perturbations to zero
     phi = 0._dp
     
  elseif (cmb_like .and. FLRW_perturb) then

     ! define filenames to read from, based on resolution
     ! !! assumes uniform grid !!
     write(deltafile,'(a,i3.3,a)')'ics_files/ics_files',res,'/delta.dat'
     write(vel1file,'(a,i3.3,a)')'ics_files/ics_files',res,'/vel1.dat'
     write(vel2file,'(a,i3.3,a)')'ics_files/ics_files',res,'/vel2.dat'
     write(vel3file,'(a,i3.3,a)')'ics_files/ics_files',res,'/vel3.dat'
     write(phifile,'(a,i3.3,a)')'ics_files/ics_files',res,'/phi.dat'

     ! open files to read phi, delta, vel perturbs from files
     open(newunit=dr_unit,file=deltafile,status='old')  ! delta rho file
     open(newunit=dv_unit1,file=vel1file,status='old')  ! delta vel file [1]
     open(newunit=dv_unit2,file=vel2file,status='old')  ! delta vel file [2]
     open(newunit=dv_unit3,file=vel3file,status='old')  ! delta vel file [3] 
     open(newunit=p_unit,file=phifile,status='old')   ! phi file

     do k = 1, cctk_gsh(3)
        do j = 1, cctk_gsh(2) ! loop over ROW no. (i is COLUMN no.) (i,j,k) --> (column, row, z)
           ! read perturbations from CMB-like generated data files
           read(dr_unit,*) delta_gs(:,j,k)
           read(dv_unit1,*) delta_vel_gs(:,j,k,1)
           read(dv_unit2,*) delta_vel_gs(:,j,k,2)
           read(dv_unit3,*) delta_vel_gs(:,j,k,3)
           read(p_unit,*), phi_gs(:,j,k)
        enddo
     enddo

     ! find indices of lower and upper bounds of location of local grid within global grid
     il = cctk_lbnd(1)
     jl = cctk_lbnd(2)
     kl = cctk_lbnd(3)
     
     iu = cctk_ubnd(1)
     ju = cctk_ubnd(2)
     ku = cctk_ubnd(3)

     ! put partial global grid into local grid
     delta = delta_gs(il:iu, jl:ju, kl:ku)
     phi = phi_gs(il:iu, jl:ju, kl:ku)
     delta_vel = delta_vel_gs(il:iu, jl:ju, kl:ku, :)

  endif

  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           if (single_mode .and. FLRW_perturb) then
              ! set single mode perturbation in phi and d(phi)/dx^i -- used for delta and delta_vel
              if (perturb_x) then
                 f = perturb_phi * sin(kx * x(i,j,k) - phi_offset)
                 df1 = perturb_phi * kx * cos(kx * x(i,j,k) - phi_offset)
                 delta_vel(i,j,k,1) = perturb_v0 * df1
              elseif (perturb_y) then
                 f = perturb_phi * sin(ky * y(i,j,k) - phi_offset)
                 df2 = perturb_phi * ky * cos(ky * y(i,j,k) - phi_offset) ! for all these cases df only has one component
                 delta_vel(i,j,k,2) = perturb_v0 * df2
              elseif (perturb_z) then
                 f = perturb_phi * sin(kz * z(i,j,k) - phi_offset)
                 df3 = perturb_phi * kz * cos(kz * z(i,j,k) - phi_offset)
                 delta_vel(i,j,k,3) = perturb_v0 * df3
              elseif (perturb_all) then
                 f = perturb_phi * (sin(kx * x(i,j,k) - phi_offset) + sin(ky * y(i,j,k) - phi_offset) + sin(kz * z(i,j,k) - phi_offset))
                 df1 = perturb_phi * kx * cos(kx * x(i,j,k) - phi_offset)
                 df2 = perturb_phi * ky * cos(ky * y(i,j,k) - phi_offset) ! df is now a vector, this is only used in velocity IC
                 df3 = perturb_phi * kz * cos(kz * z(i,j,k) - phi_offset)
                 delta_vel(i,j,k,1) = perturb_v0 * df1
                 delta_vel(i,j,k,2) = perturb_v0 * df2
                 delta_vel(i,j,k,3) = perturb_v0 * df3
              endif
              phi(i,j,k) = f
              phidot = 0._dp     !  d(phi)/dt
              delta(i,j,k) = perturb_rho0 * f              
           endif
       
           !! set up metric, extrinsic curvature, lapse and shift
           if (data) then

              if (lapse) then
                 !if (FLRW_perturb_metric) then
                 if (FLRW_perturb) then
                    alp(i,j,k) = sqrt(1._dp + 2._dp * phi(i,j,k))
                 else
                    alp(i,j,k) = FLRW_lapse_value
                 endif
              endif
                 
              !if (FLRW_perturb_metric) then        
              if (FLRW_perturb) then
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

              ! May also need the derivative of the lapse -- evolution of this is specified in ADMBase.
              if (dtlapse) then
                 dtalp(i,j,k) = 0._dp
              end if              

              if (shift) then
                 betax(i,j,k) = 0._dp
                 betay(i,j,k) = 0._dp
                 betaz(i,j,k) = 0._dp
              end if
         
              ! set up  matter variables
              if (hydro) then
                 ! initialise matter to homogeneous values
                 rho(i,j,k) = rho0
                 !           press(i,j,k) = FLRW_K * rho0 ** FLRW_gamma   !! previously. removing these FLRW keywords
                 press(i,j,k) = 0._dp
                 eps(i,j,k) = 0._dp
                 vel(i,j,k,:) = 0._dp

!                 if (FLRW_perturb_density) then
                 if (FLRW_perturb) then

                    vel(i,j,k,:) = delta_vel(i,j,k,:)
                    rho(i,j,k) = rho(i,j,k) * ( 1._dp + delta(i,j,k) )

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
