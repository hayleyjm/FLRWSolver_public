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
  real :: a0, kvalue, asq, adot, rho0, r_gauss, r0, perturb_rho0, box_length, phi, rad, kx
  real, parameter :: pi = 3.14159265358979323846264338327
  real :: P, Q, W, offset_x, offset_y, offset_z
  logical   :: lapse, dtlapse, shift, data, hydro
  lapse = CCTK_EQUALS (initial_lapse, "flrw")
  dtlapse = CCTK_EQUALS (initial_dtlapse, "flrw")
  shift = CCTK_EQUALS (initial_shift, "flrw")
  data  = CCTK_EQUALS (initial_data,  "flrw")
  hydro = CCTK_EQUALS (initial_hydro, "flrw")

  rho0 = FLRW_init_rho
  perturb_rho0 = FLRW_pert_amplitude*rho0
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

  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)

	 !! set phi depending on perturbation

	 if (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then

	    if (CCTK_EQUALS (FLRW_phi_type, "General")) then
	       phi = perturb_rho0*sin(kx*x(i,j,k))
	    elseif (CCTK_EQUALS (FLRW_phi_type, "Poisson")) then	    
	       phi = -4.*pi*perturb_rho0*sin(kx*x(i,j,k))/kx**2
	    endif

	 elseif (CCTK_EQUALS (FLRW_perturb_type, "Tophat")) then

	 	rad = sqrt((x(i,j,k) - offset_x)**2 + (y(i,j,k) - offset_y)**2 + (z(i,j,k) - offset_z)**2)

		if (rad <= r0) then
         	   phi = 2.*pi* (rho0 + perturb_rho0)*(rad**2 - 3.*r0**2)/3.
		else
		   phi = -4.*pi*r0**2*rho0/3.
		endif
	endif
        
	 !! set up metric, extrinsic curvature, lapse and shift

	 if (data) then
	    
	    if (FLRW_perturb_metric) then
              gxx(i,j,k) = asq*(1.0 - 2.0*phi)
              gxy(i,j,k) = 0.0
              gxz(i,j,k) = 0.0
              gyy(i,j,k) = asq*(1.0 - 2.0*phi)
              gyz(i,j,k) = 0.0
              gzz(i,j,k) = asq*(1.0 - 2.0*phi)
            
              kxx(i,j,k) = kvalue*(1.0 - 2.0*phi)/sqrt(1.0 + 2.0*phi)
              kxy(i,j,k) = 0.0
              kxz(i,j,k) = 0.0
              kyy(i,j,k) = kvalue*(1.0 - 2.0*phi)/sqrt(1.0 + 2.0*phi)
              kyz(i,j,k) = 0.0
              kzz(i,j,k) = kvalue*(1.0 - 2.0*phi)/sqrt(1.0 + 2.0*phi)
            
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
	 end if

         if (lapse) then
	    if (FLRW_perturb_metric) then
               alp(i,j,k) = sqrt(1.0 + 2.0*phi)
   	    else
	       alp(i,j,k) = 1.0
            endif
	 end if

         ! May also need the derivative of the lapse -- this is specified in ADMBase (somehow).
         if (dtlapse) then
            dtalp(i,j,k)= 0.0
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
            vel(i,j,k,1) = 0.0
            vel(i,j,k,2) = 0.0
            vel(i,j,k,3) = 0.0
         end if
         
	 !! set density depending on perturbation

	if (FLRW_perturb_density) then

	   if (CCTK_EQUALS (FLRW_perturb_type, "Gaussian")) then	

	      	rad = sqrt((x(i,j,k) - offset_x)**2 + (y(i,j,k) - offset_y)**2 + (z(i,j,k) - offset_z)**2)
		rho(i,j,k) = rho(i,j,k) + perturb_rho0*exp(-rad**2/r0**2)

	   elseif (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then
    	
	   	if (CCTK_EQUALS (FLRW_phi_type, "General")) then   

		   	P = perturb_rho0**2*kx**2
	    	   	Q = perturb_rho0*kx**2
            	   	W = sqrt(2.*pi*rho0/3.)

	    	   	rho(i,j,k) = ((6.*(-2.*W + 4.*perturb_rho0*W*sin(kx*x(i,j,k)))**2) / ((-1. + 2.*perturb_rho0*sin(kx*x(i,j,k)))**2 * (1. + 2.*perturb_rho0*sin(kx*x(i,j,k)))) + (2.*(-3.*P*cos(kx*x(i,j,k))**2 + 2.*Q*sin(kx*x(i,j,k)) - 4.*P*sin(kx*x(i,j,k))**2)) / ((-1 + 2.*perturb_rho0*sin(kx*x(i,j,k)))**3)) / (16*pi)

		elseif (CCTK_EQUALS (FLRW_phi_type, "Poisson")) then

	            rho(i,j,k) = rho(i,j,k) + perturb_rho0*sin(kx*x(i,j,k))

		endif

	    elseif (CCTK_EQUALS (FLRW_perturb_type, "Tophat")) then

	    	if (rad < r0) then
		   rho(i,j,k) = rho(i,j,k) + perturb_rho0
		else
		   rho(i,j,k) = rho(i,j,k)
		endif   

	    endif

	    ! COPY TEST DENSITY HERE - SO AS TO NOT DISTURB OTHER DENSITIES
	    !! current test: density/4pi

	    if (FLRW_test) then

	       rho(i,j,k) = ((6.*(-2.*W + 4.*perturb_rho0*W*sin(kx*x(i,j,k)))**2) / ((-1. + 2.*perturb_rho0*sin(kx*x(i,j,k)))*\
*2 * (1. + 2.*perturb_rho0*sin(kx*x(i,j,k)))) + (2.*(-3.*P*cos(kx*x(i,j,k))**2 + 2.*Q*sin(kx*x(i,j,k)) - 4.*P*sin(kx*x(i,j,k))**2)) / ((-1 + 2.*perturb_rho0*sin(kx*x(i,j,k)))**3)/(16*pi) - rho0)/(4.0*pi) + rho0

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
