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
  real(dp) :: a0, kvalue, asq, adot, rho0, r_gauss, r0, perturb_rho0, box_length_x, box_length_y, box_length_z, phi, rad, kx, perturb_phi, perturb_rho0_rel, delphi, delsqphi, phidot
  real(dp), parameter :: pi = 3.14159265358979323846264338327
  real(dp) :: P, Q, W, offset_x, offset_y, offset_z, lapse_value, phi_offset, amp, perturb_v0, ky, kz, f, df, H0, z0, df1, df2, df3
  logical   :: lapse, dtlapse, shift, data, hydro, perturb_x, perturb_y, perturb_z, perturb_all
  lapse = CCTK_EQUALS (initial_lapse, "flrw")
  dtlapse = CCTK_EQUALS (initial_dtlapse, "flrw")
  shift = CCTK_EQUALS (initial_shift, "flrw")
  data  = CCTK_EQUALS (initial_data,  "flrw")
  hydro = CCTK_EQUALS (initial_hydro, "flrw")
  perturb_x = CCTK_EQUALS (FLRW_perturb_direction, "x")
  perturb_y = CCTK_EQUALS (FLRW_perturb_direction, "y")
  perturb_z = CCTK_EQUALS (FLRW_perturb_direction, "z")
  perturb_all = CCTK_EQUALS (FLRW_perturb_direction, "all")

  z0 = FLRW_initial_redshift

! testing initial conds for a, rho, calculated from "real" length scale etc
!  if (FLRW_test_ics) then
!     H0 = 5._dp / 30._dp ! hubble parameter in code units
!     a0 = 1._dp / (1._dp + z0) ! scale factor found from initial redshift of z=127 (as in Milennium)
!     rho0 = 3._dp * H0**2 / (8._dp * pi * a0**3) ! rho = rhoc * a**-3
!  else ! regular IC's
  rho0 = FLRW_init_rho
  a0 =1._dp / (1._dp + z0)
!  endif

  amp = FLRW_phi_pert_amplitude
  perturb_phi = amp * rho0
  r0 = FLRW_radius
  phi_offset = FLRW_phi_phase_offset
  asq = a0*a0
  adot = sqrt((8.0_dp * pi)*rho0*asq/3.0_dp)
  kvalue = -adot * a0
  box_length_x = 2.0_dp * FLRW_xmax
  box_length_y = 2.0_dp * FLRW_xmax !! these are all equal for now, will need to introduce FLRW_y/zmax if we choose uneqal box
  box_length_z = 2.0_dp * FLRW_xmax
  offset_x = FLRW_offset_x
  offset_y = FLRW_offset_y
  offset_z = FLRW_offset_z
  kx = 2.0_dp*pi/box_length_x
  ky = 2.0_dp*pi/box_length_y
  kz = 2.0_dp*pi/box_length_z
  lapse_value = FLRW_lapse_value	!! only needed for FLRW

!############################################### TEST PRINT STATEMENT
  print*,'SETTING UP FLRW - YOUR PRINT STATEMENT WORKS!'
! ##################################################################

! set density, velocity amplitudes
  if (CCTK_EQUALS (FLRW_phi_solution, "Constant")) then !constant(phi) mode
      
     perturb_rho0 = - (kx**2 / (4.0_dp * pi * rho0 * asq) + 2.0_dp)
!     perturb_rho0 = - (kx**2 / (4.0_dp * pi * rho0 * asq))
     perturb_v0 = - 1._dp / (asq * sqrt(6._dp * pi * rho0))

  elseif (CCTK_EQUALS (FLRW_phi_solution, "Decaying")) then !decaying mode
  
     perturb_rho0 = (3.0_dp * kx**2 / (20.0_dp * pi * rho0 * asq) - 9.0_dp / 5.0_dp)
     perturb_v0 = - 3._dp * sqrt(3._dp / 8._dp * pi * rho0) / (5._dp * asq)
 
  endif

  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)

         ! set f(or g) for perturbations
         if (perturb_x) then
            f = perturb_phi * sin(kx * x(i,j,k) - phi_offset)
            df1 = perturb_phi * kx * cos(kx * x(i,j,k) - phi_offset)
         elseif (perturb_y) then
            f = perturb_phi * sin(ky* y(i,j,k) - phi_offset)
            df2 = perturb_phi * ky * cos(ky * y(i,j,k) - phi_offset) ! for all these cases df only has one component
         elseif (perturb_z) then
            f = perturb_phi * sin(kz* z(i,j,k) - phi_offset)
            df3 = perturb_phi * kz * cos(kz * z(i,j,k) - phi_offset)
         elseif (perturb_all) then
            f = perturb_phi * (sin(kx* x(i,j,k) - phi_offset) + sin(ky* y(i,j,k) - phi_offset) + sin(kz* z(i,j,k) - phi_offset))
            df1 = perturb_phi * kx * cos(kx * x(i,j,k) - phi_offset)
            df2 = perturb_phi * ky * cos(ky * y(i,j,k) - phi_offset) ! df is now a vector, this is only used in velocity IC
            df3 = perturb_phi * kz * cos(kz * z(i,j,k) - phi_offset)
         endif

	 !! set phi depending on perturbation type
	 if (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then

               if (CCTK_EQUALS (FLRW_phi_solution, "Constant")) then !constant(phi) solution

                  phi = f
!                  delphi = df ! don't need these right now, use df anyway in velocity
 !                 delsqphi = -kx**2 * phi
                  phidot = 0._dp
               
               elseif (CCTK_EQUALS (FLRW_phi_solution, "Decaying")) then !decaying solution
               
                  phi = -3._dp / 5._dp * f
!                  delphi = -3._dp / 5._dp * df ! don't use these, use df instead
!                  delsqphi = -kx**2 * phi
                  phidot = sqrt(6._dp * pi * rho0) * f
               
               endif
	endif
        
	 !! set up metric, extrinsic curvature, lapse and shift
	 if (data) then

	    if (lapse) then
               if (FLRW_perturb_metric) then

               	  alp(i,j,k) = sqrt(1.0_dp + 2.0_dp * phi)

               else

		  alp(i,j,k) = lapse_value

               endif
            end if

	   
	    if (FLRW_perturb_metric) then

             if (FLRW_negative_metric) then
                gxx(i,j,k) = -asq*(1.0_dp - 2.0_dp*phi)
                gxy(i,j,k) = 0.0
                gxz(i,j,k) = 0.0
                gyy(i,j,k) = -asq*(1.0_dp - 2.0_dp*phi)
                gyz(i,j,k) = 0.0
                gzz(i,j,k) = -asq*(1.0_dp - 2.0_dp*phi)
             else
                gxx(i,j,k) = asq*(1.0_dp - 2.0_dp*phi)
                gxy(i,j,k) = 0.0
                gxz(i,j,k) = 0.0
                gyy(i,j,k) = asq*(1.0_dp - 2.0_dp*phi)
                gyz(i,j,k) = 0.0
                gzz(i,j,k) = asq*(1.0_dp - 2.0_dp*phi)
             endif

              if (FLRW_negative_curv) then
                 kxx(i,j,k) = -kvalue*(1.0_dp - 2.0_dp * phi - phidot * a0 / adot) / alp(i,j,k)
                 kxy(i,j,k) = 0.0
                 kxz(i,j,k) = 0.0
                 kyy(i,j,k) = -kvalue*(1.0_dp - 2.0_dp * phi - phidot * a0 / adot) / alp(i,j,k)
                 kyz(i,j,k) = 0.0
                 kzz(i,j,k) = -kvalue*(1.0_dp - 2.0_dp * phi - phidot * a0 / adot) / alp(i,j,k)
              else
                 kxx(i,j,k) = kvalue*(1.0_dp - 2.0_dp * phi - phidot * a0 / adot) / alp(i,j,k)
                 kxy(i,j,k) = 0.0
                 kxz(i,j,k) = 0.0
                 kyy(i,j,k) = kvalue*(1.0_dp - 2.0_dp * phi - phidot * a0 / adot) / alp(i,j,k)
                 kyz(i,j,k) = 0.0
                 kzz(i,j,k) = kvalue*(1.0_dp - 2.0_dp * phi - phidot * a0 / adot) / alp(i,j,k)
              endif

	    else

              if (FLRW_negative_metric) then
                 gxx(i,j,k) = -asq
                 gxy(i,j,k) = 0.0
                 gxz(i,j,k) = 0.0
                 gyy(i,j,k) = -asq
                 gyz(i,j,k) = 0.0
                 gzz(i,j,k) = -asq
              else
                 gxx(i,j,k) = asq
                 gxy(i,j,k) = 0.0
                 gxz(i,j,k) = 0.0
                 gyy(i,j,k) = asq
                 gyz(i,j,k) = 0.0
                 gzz(i,j,k) = asq
              endif

              kxx(i,j,k) = kvalue
              kxy(i,j,k) = 0.0
              kxz(i,j,k) = 0.0
              kyy(i,j,k) = kvalue
              kyz(i,j,k) = 0.0
              kzz(i,j,k) = kvalue
	    endif
	 endif

         ! May also need the derivative of the lapse -- this is specified in ADMBase (somehow).
         if (dtlapse) then
	    dtalp(i,j,k) = 0.0
	 end if


         if (shift) then
            betax(i,j,k) = 0.0
            betay(i,j,k) = 0.0
            betaz(i,j,k) = 0.0
         end if
         
         ! set up  matter variables
         if (hydro) then
            rho(i,j,k) = rho0
            press(i,j,k) = FLRW_K * rho0 ** FLRW_gamma
            eps(i,j,k) = 0.0
            if (FLRW_perturb_density) then

               if (perturb_x) then
                  vel(i,j,k,1) = perturb_v0 * df1
                  vel(i,j,k,2) = 0.0
                  vel(i,j,k,3) = 0.0
               elseif (perturb_y) then
                  vel(i,j,k,1) = 0.0
                  vel(i,j,k,2) = perturb_v0 * df2
                  vel(i,j,k,3) = 0.0
               elseif (perturb_z) then
                  vel(i,j,k,1) = 0.0
                  vel(i,j,k,2) = 0.0
                  vel(i,j,k,3) = perturb_v0 * df3
               elseif (perturb_all) then
                  vel(i,j,k,1) = perturb_v0 * df1
                  vel(i,j,k,2) = perturb_v0 * df2
                  vel(i,j,k,3) = perturb_v0 * df3
               endif

            else
               vel(i,j,k,1) = 0.0
               vel(i,j,k,2) = 0.0
               vel(i,j,k,3) = 0.0
            end if
         end if
         
	 !! set density depending on perturbation

	if (FLRW_perturb_density) then

	   if (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then
               
              rho(i,j,k) = rho(i,j,k) * ( 1._dp + perturb_rho0 * f)
           
           endif

	endif

        end do
     end do
  end do

  if (CCTK_EQUALS (metric_type, "physical")) then
     ! do nothing
  else
     call CCTK_WARN (0, "Unknown value of ADMBase::metric_type -- FLRW only set-up for metric_type = physical")
  end if

end subroutine FLRW_InitialData
