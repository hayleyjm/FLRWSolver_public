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
  real :: a0, kvalue, asq, adot, rho0, r_gauss, r0, perturb_rho0, box_length, phi, tophat_rad, rad, kx
  real, parameter :: pi = 3.14159265358979323846264338327
  logical   :: lapse, dtlapse, shift, data, hydro
  lapse = CCTK_EQUALS (initial_lapse, "flrw")
  dtlapse = CCTK_EQUALS (initial_dtlapse, "flrw")
  shift = CCTK_EQUALS (initial_shift, "flrw")
  data  = CCTK_EQUALS (initial_data,  "flrw")
  hydro = CCTK_EQUALS (initial_hydro, "flrw")

  rho0 = FLRW_init_rho
  perturb_rho0 = FLRW_pert_amplitude
  r0 = Gaussian_r0
  a0 = 1.0
  asq = a0*a0
  adot = sqrt((8.0 * pi)*rho0*asq/3.0)
  kvalue = -adot * a0
  box_length = 480
  tophat_rad = 150
  kx = 2.0*pi/box_length

  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)

         rad = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)

	 !! set phi depending on perturbation

	 if (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then
	    phi = -perturb_rho0*box_length**2*sin(kx*x(i,j,k))/pi
	 elseif (CCTK_EQUALS (FLRW_perturb_type, "Tophat")) then
	 	if (rad <= tophat_rad) then
         	   phi = -0.25*perturb_rho0*(tophat_rad**2 - rad**2)
		else
		   phi = 0.
		endif
	endif
        
	 !! set up metric, extrinsic curvature, lapse and shift

	 if (data) then
            gxx(i,j,k) = asq*(1.0 - 2.0*phi)
            gxy(i,j,k) = 0.0
            gxz(i,j,k) = 0.0
            gyy(i,j,k) = asq*(1.0 - 2.0*phi)
            gyz(i,j,k) = 0.0
            gzz(i,j,k) = asq*(1.0 - 2.0*phi)
            
            kxx(i,j,k) = kvalue*(1.0 - 2.0*phi)/sqrt(1.0 + 2.0*phi)
            kxy(i,j,k) = 0.0
            kxz(i,j,k) = 0.0
            kyy(i,j,k) = kvalue*(1.0 - 3.0*phi)/sqrt(1.0 + 2.0*phi)
            kyz(i,j,k) = 0.0
            kzz(i,j,k) = kvalue*(1.0 - 3.0*phi)/sqrt(1.0 + 2.0*phi)
         end if

         if (lapse) then
            alp(i,j,k) = sqrt(1.0 + 2.0*phi)
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
         
	if (FLRW_perturb) then
	    if (CCTK_EQUALS (FLRW_perturb_type, "Gaussian")) then
		r_gauss = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
		rho(i,j,k) = rho(i,j,k) + perturb_rho0*exp(-r_gauss**2/r0**2)
	    elseif (CCTK_EQUALS (FLRW_perturb_type, "Sine")) then
	        rho(i,j,k) = rho(i,j,k) + perturb_rho0*sin(kx*x(i,j,k))
	    elseif (CCTK_EQUALS (FLRW_perturb_type, "Tophat")) then
	    	if (rad < tophat_rad) then
		   rho(i,j,k) = rho(i,j,k) + perturb_rho0
		else
		   rho(i,j,k) = rho(i,j,k)
		endif   
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
