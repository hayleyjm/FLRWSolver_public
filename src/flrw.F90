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
  real :: a0, kvalue, asq, adot, rho0, r_gauss, r0, perturb_rho0, box_length, phi, rad, kx, perturb_phi, perturb_rho0_rel, delphi, delsqphi, phidot
  real, parameter :: pi = 3.14159265358979323846264338327
  real :: P, Q, W, offset_x, offset_y, offset_z, lapse_value
  logical   :: lapse, dtlapse, shift, data, hydro
  lapse = CCTK_EQUALS (initial_lapse, "flrw")
  dtlapse = CCTK_EQUALS (initial_dtlapse, "flrw")
  shift = CCTK_EQUALS (initial_shift, "flrw")
  data  = CCTK_EQUALS (initial_data,  "flrw")
  hydro = CCTK_EQUALS (initial_hydro, "flrw")

  rho0 = FLRW_init_rho
  perturb_phi = FLRW_phi_pert_amplitude * rho0
  r0 = FLRW_radius
  a0 = 1.0
  asq = a0*a0
  adot = sqrt((8.0 * pi)*rho0*asq/3.0)
  kvalue = -adot * a0
  box_length = 2.0 * FLRW_xmax
  offset_x = FLRW_offset_x
  offset_y = FLRW_offset_y
  offset_z = FLRW_offset_z
  kx = 2.0*pi/box_length
  perturb_rho0 = -kx**2 * perturb_phi / (4.0 * pi)
  perturb_rho0_rel = - perturb_phi * (kx**2 + 3.0 * adot**2) / (4.0 * pi * asq)
  lapse_value = FLRW_lapse_value	!! only needed for FLRW


  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)

	 !! set phi depending on perturbation

	 if (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then

	    	!! phi same for both newtonian & relativistic theory cases now

	       phi = perturb_phi * sin(kx * x(i,j,k))
               delphi = kx * perturb_phi * cos(kx * x(i,j,k))
               delsqphi = -kx**2 * phi
               phidot = 0.0

	 elseif (CCTK_EQUALS (FLRW_perturb_type, "Tophat")) then

	 	rad = sqrt((x(i,j,k) - offset_x)**2 + (y(i,j,k) - offset_y)**2 + (z(i,j,k) - offset_z)**2)

		if (rad <= r0) then
         	   phi = 2.0 * pi * (rho0 + perturb_rho0) * (rad**2 - 3.0 * r0**2) / 3.0
		else
		   phi = -4.0 * pi * r0**2 * rho0 / 3.0
		endif
	endif
        
	 !! set up metric, extrinsic curvature, lapse and shift

	 if (data) then

	    if (lapse) then
               if (FLRW_perturb_metric) then

               	  alp(i,j,k) = sqrt(1.0 + 2.0 * phi)

               else

		  alp(i,j,k) = lapse_value

               endif
            end if

	   
	    if (FLRW_perturb_metric) then
              gxx(i,j,k) = asq*(1.0 - 2.0*phi)
              gxy(i,j,k) = 0.0
              gxz(i,j,k) = 0.0
              gyy(i,j,k) = asq*(1.0 - 2.0*phi)
              gyz(i,j,k) = 0.0
              gzz(i,j,k) = asq*(1.0 - 2.0*phi)
            
              kxx(i,j,k) = kvalue*(1.0 - 2.0 * phi)/alp(i,j,k)
              kxy(i,j,k) = 0.0
              kxz(i,j,k) = 0.0
              kyy(i,j,k) = kvalue*(1.0 - 2.0 * phi)/alp(i,j,k)
              kyz(i,j,k) = 0.0
              kzz(i,j,k) = kvalue*(1.0 - 2.0 * phi)/alp(i,j,k)

	    else
	      gxx(i,j,k) = asq
              gxy(i,j,k) = 0.0
              gxz(i,j,k) = 0.0
              gyy(i,j,k) = asq
              gyz(i,j,k) = 0.0
              gzz(i,j,k) = asq

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
            vel(i,j,k,1) = -delphi * adot / (asq*a0 * 4.0 * pi * rho0)
            vel(i,j,k,2) = 0.0
            vel(i,j,k,3) = 0.0
         end if
         
	 !! set density depending on perturbation

	if (FLRW_perturb_density) then

	   if (CCTK_EQUALS (FLRW_perturb_type, "Gaussian")) then	

	      	rad = sqrt((x(i,j,k) - offset_x)**2 + (y(i,j,k) - offset_y)**2 + (z(i,j,k) - offset_z)**2)
		rho(i,j,k) = rho(i,j,k) + perturb_rho0*exp(-rad**2/r0**2)

	   elseif (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then
   
		if (CCTK_EQUALS (FLRW_density_type, "Hamiltonian")) then		!!set general rho based on phi

		   	rho(i,j,k) = ((2.0 * (-3.0 * delphi**2 - 2.0 * delsqphi + 4.0 * phi * delsqphi)) / (asq * (-1.0 + 2.0 * phi)**3) + (6.0 * (-adot + 2.0 * phi * adot + a0 * phidot)**2) / (asq * (-1.0 + 2.0 * phi)**2 * alp(i,j,k)**2)) / (16.0 * pi)


		elseif (CCTK_EQUALS (FLRW_perturb_theory_type, "Relativistic") .and. CCTK_EQUALS (FLRW_density_type, "Poisson")) then

	               rho(i,j,k) = rho(i,j,k) + perturb_rho0_rel * sin(kx * x(i,j,k))

                elseif (CCTK_EQUALS (FLRW_perturb_theory_type, "Newtonian") .and. CCTK_EQUALS (FLRW_density_type, "Poisson")) then

                       rho(i,j,k) = rho(i,j,k) + perturb_rho0 * sin(kx * x(i,j,k))
		
		endif

	    elseif (CCTK_EQUALS (FLRW_perturb_type, "Tophat")) then

	    	if (rad < r0) then
		   rho(i,j,k) = rho(i,j,k) + perturb_rho0
		else
		   rho(i,j,k) = rho(i,j,k)
		endif   

	    endif




	    ! COPY TEST DENSITY HERE - SO AS TO NOT DISTURB OTHER DENSITIES
	    !! current tests: density/4pi and different gauge alpha (which will make no difference)

	    if (FLRW_test) then

!!	       rho(i,j,k) = ((6.*(-2.*W + 4.*perturb_rho0*W*sin(kx*x(i,j,k)))**2) / ((-1. + 2.*perturb_rho0*sin(kx*x(i,j,k)))*\
*2 * (1. + 2.*perturb_rho0*sin(kx*x(i,j,k)))) + (2.*(-3.*P*cos(kx*x(i,j,k))**2 + 2.*Q*sin(kx*x(i,j,k)) - 4.*P*sin(kx*x(i,j,k))**2)) / ((-1 + 2.*perturb_rho0*sin(kx*x(i,j,k)))**3)/(16*pi) - rho0)/(4.0*pi) + rho0

     	       alp(i,j,k) = a0 * sqrt(1 + 2.0 * phi)

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




!! OLD HAMILTONIAN DENSITIES -- BEFORE I PUT IN GENERAL
!! DO NOT DELETE JUST YET - UNTIL I AM SURE THE GENERAL ONE WORKS.

!!!             if (CCTK_EQUALS (FLRW_perturb_theory_type, "Relativistic") .and. CCTK_EQUALS (FLRW_d\
ensity_type, "Hamiltonian")) then

!!!                     print*,' SETTING RELATIVISTIC HAMILTONIAN DENSITY'

!!!                     P = perturb_phi**2 * kx**2 / asq
!!!                     Q = perturb_phi * kx**2 / a0
!!!                     W = perturb_phi / a0

!!!                     rho(i,j,k) = ((2.0 * (-3.0 * P * cos(kx*x(i,j,k))**2 + 2.0 * Q * sin(kx*x(i,\
j,k)) - 4.0 * P * sin(kx*x(i,j,k))**2)) / (asq * (-1.0 + 2.0 * W * sin(kx*x(i,j,k)))**3) + (6.0 * (-\
adot + W * adot * sin(kx*x(i,j,k)))**2) / (asq * (-1.0 + 2.0 * W * sin(kx*x(i,j,k)))**2 * alp(i,j,k)\
**2)) / (16.0 * pi)





!!!             elseif (CCTK_EQUALS (FLRW_perturb_theory_type, "Newtonian") .and. CCTK_EQUALS (FLRW_\
density_type, "Hamiltonian")) then

!!!                    print*,'SETTING NEWTONIAN HAMILTONIAN DENSITY'

!!!                     P = perturb_phi**2 * kx**2
!!!                        Q = perturb_phi * kx**2
!!!                        W = sqrt(2.0 * pi * rho0 / 3.0)

!!!                rho(i,j,k) = ((6.0 * (-2.0 * W + 4.0 * perturb_phi * W * sin(kx * x(i,j,k)))**2) \
/ ((-1.0 + 2.0 * perturb_phi * sin(kx*x(i,j,k)))**2 * (1.0 + 2.0 * perturb_phi * sin(kx * x(i,j,k)))\
) + (2.0 * (-3.0 * P * cos(kx * x(i,j,k))**2 + 2.0 * Q * sin(kx * x(i,j,k)) - 4.0 * P * sin(kx * x(i\
,j,k))**2)) / ((-1.0 + 2.0 * perturb_phi * sin(kx * x(i,j,k)))**3)) / (16.0 * pi)













end subroutine FLRW_InitialData
