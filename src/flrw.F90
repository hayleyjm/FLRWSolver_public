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
  real(dp) :: a0, kvalue, asq, adot, rho0, r_gauss, r0, perturb_rho0, box_length, phi, rad, kx, perturb_phi, perturb_rho0_rel, delphi, delsqphi, phidot
  real(dp), parameter :: pi = 3.14159265358979323846264338327
  real(dp) :: P, Q, W, offset_x, offset_y, offset_z, lapse_value, phi_offset, amp, perturb_v0
  logical   :: lapse, dtlapse, shift, data, hydro
  lapse = CCTK_EQUALS (initial_lapse, "flrw")
  dtlapse = CCTK_EQUALS (initial_dtlapse, "flrw")
  shift = CCTK_EQUALS (initial_shift, "flrw")
  data  = CCTK_EQUALS (initial_data,  "flrw")
  hydro = CCTK_EQUALS (initial_hydro, "flrw")

  rho0 = FLRW_init_rho
  amp = FLRW_phi_pert_amplitude
  perturb_phi = amp * rho0
  r0 = FLRW_radius
  phi_offset = FLRW_phi_phase_offset
  a0 = 1.0_dp
  asq = a0*a0
  adot = sqrt((8.0_dp * pi)*rho0*asq/3.0_dp)
  kvalue = -adot * a0
  box_length = 2.0_dp * FLRW_xmax
  offset_x = FLRW_offset_x
  offset_y = FLRW_offset_y
  offset_z = FLRW_offset_z
  kx = 2.0_dp*pi/box_length
 !! perturb_rho0 = -kx**2 * perturb_phi / (4.0_dp * pi)
  lapse_value = FLRW_lapse_value	!! only needed for FLRW

  if (CCTK_EQUALS (FLRW_phi_solution, "Constant")) then !constant(phi) mode
      
  !!   perturb_rho0_rel = - perturb_phi * (kx**2 + 3.0_dp * adot**2) / (4.0_dp * pi * asq)
     perturb_rho0 = - rho0 * (kx**2 * amp / (4.0_dp * pi * asq) + 2.0_dp * amp * rho0)
     perturb_v0 = - 1._dp / (asq * sqrt(6._dp * pi * rho0)) !! original velocity
!!     perturb_v0 = -10._dp * sqrt(6._dp) / (24._dp * asq * sqrt(pi * rho0))  !! velocity from mom constraint with covD(K_j^i)

  elseif (CCTK_EQUALS (FLRW_phi_solution, "Decaying")) then !decaying mode
! decaying soln needs to be checked before running!!!
  
!!     perturb_rho0_rel = perturb_phi * (3.0_dp * kx**2 / (20._dp * pi) - 9._dp * rho0 / 5._dp)
     perturb_rho0 = perturb_phi * (3.0_dp * kx**2 / (20.0_dp * pi * rho0 * a0**2) - 9.0_dp / 5.0_dp)
 
  endif

  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)

	 !! set phi depending on perturbation type

	 if (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then

               if (CCTK_EQUALS (FLRW_phi_solution, "Constant")) then !constant(phi) solution
               
                  phi = perturb_phi * sin(kx * x(i,j,k) - phi_offset)
                  delphi = kx * perturb_phi * cos(kx * x(i,j,k) - phi_offset)
                  delsqphi = -kx**2 * phi
                  phidot = 0._dp
               
               elseif (CCTK_EQUALS (FLRW_phi_solution, "Decaying")) then !decaying solution
               
                  phi = -3._dp / 5._dp * perturb_phi * sin(kx * x(i,j,k) - phi_offset)
                  delphi = -3._dp / 5._dp * kx * perturb_phi * cos(kx*x(i,j,k) - phi_offset)
                  delsqphi = -kx**2 * phi
                  phidot = sqrt(6._dp * pi * rho0) * perturb_phi * sin(kx * x(i,j,k) - phi_offset)
               
               endif
               
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

               	  alp(i,j,k) = sqrt(1.0_dp + 2.0_dp * phi)

               else

		  alp(i,j,k) = lapse_value

               endif
            end if

	   
	    if (FLRW_perturb_metric) then
              gxx(i,j,k) = asq*(1.0_dp - 2.0_dp*phi)
              gxy(i,j,k) = 0.0
              gxz(i,j,k) = 0.0
              gyy(i,j,k) = asq*(1.0_dp - 2.0_dp*phi)
              gyz(i,j,k) = 0.0
              gzz(i,j,k) = asq*(1.0_dp - 2.0_dp*phi)
 
              kxx(i,j,k) = kvalue*(1.0_dp - 2.0_dp * phi - phidot * a0 / adot) / alp(i,j,k)
              kxy(i,j,k) = 0.0
              kxz(i,j,k) = 0.0
              kyy(i,j,k) = kvalue*(1.0_dp - 2.0_dp * phi - phidot * a0 / adot) / alp(i,j,k)
              kyz(i,j,k) = 0.0
              kzz(i,j,k) = kvalue*(1.0_dp - 2.0_dp * phi - phidot * a0 / adot) / alp(i,j,k)

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
            if (FLRW_perturb_density) then

               if (CCTK_EQUALS (FLRW_phi_solution, "Constant")) then ! constant(phi) mode

                  vel(i,j,k,1) = perturb_v0 * delphi
                  vel(i,j,k,2) = 0.0
                  vel(i,j,k,3) = 0.0

               elseif (CCTK_EQUALS (FLRW_phi_solution, "Decaying")) then !decaying mode
                  ! decaying soln needs to be checked before running !!!
                  vel(i,j,k,1) = - 3._dp * delphi * sqrt(3._dp / (32._dp * pi * rho0)) / (a0**2 * 5._dp)
                  vel(i,j,k,2) = 0.0
                  vel(i,j,k,3) = 0.0

               endif

            else
               vel(i,j,k,1) = 0.0
               vel(i,j,k,2) = 0.0
               vel(i,j,k,3) = 0.0
            end if
         end if
         
	 !! set density depending on perturbation

	if (FLRW_perturb_density) then

	   if (CCTK_EQUALS (FLRW_perturb_type, "Gaussian")) then	

	      	rad = sqrt((x(i,j,k) - offset_x)**2 + (y(i,j,k) - offset_y)**2 + (z(i,j,k) - offset_z)**2)
		rho(i,j,k) = rho(i,j,k) + perturb_rho0*exp(-rad**2/r0**2)

	   elseif (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then
   
                    rho(i,j,k) = rho(i,j,k) + perturb_rho0 * sin(kx * x(i,j,k) - phi_offset)
                    

!                if (CCTK_EQUALS (FLRW_density_type, "Hamiltonian")) then     !!set general rho based on phi

!		   	rho(i,j,k) = ((2.0 * (-3.0 * delphi**2 - 2.0 * delsqphi + 4.0 * phi * delsqphi)) / (asq * (-1.0 + 2.0 * phi)**3) + (6.0 * (-adot + 2.0 * phi * adot + a0 * phidot)**2) / (asq * (-1.0 + 2.0 * phi)**2 * alp(i,j,k)**2)) / (16.0 * pi)


!		elseif (CCTK_EQUALS (FLRW_perturb_theory_type, "Relativistic") .and. CCTK_EQUALS (FLRW_density_type, "Poisson")) then

!	               rho(i,j,k) = rho(i,j,k) + perturb_rho0_rel * sin(kx * x(i,j,k) - phi_offset)

 !               elseif (CCTK_EQUALS (FLRW_perturb_theory_type, "Newtonian") .and. CCTK_EQUALS (FLRW_density_type, "Poisson")) then

  !                     rho(i,j,k) = rho(i,j,k) + perturb_rho0 * sin(kx * x(i,j,k) - phi_offset)
		
!		endif

                    ! tophat doesnt work
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

end subroutine FLRW_InitialData
