! file    flrw_singlemode.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    Single-mode perturbation of FLRW initial data, longitudinal gauge, zero shift
!            note: these ICs assume linear perturbation theory is valid --> O(10^-10) violation of constraints

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine FLRW_SingleMode (CCTK_ARGUMENTS)
  USE FLRW_InitTools
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer   :: i,j,k
  logical   :: lapse,dtlapse,shift,data,hydro
  logical   :: perturb_x,perturb_y,perturb_z,perturb_all

  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
  CCTK_REAL :: lambda(3),kval,kx,ky,kz,perturb_rho0,perturb_v0
  CCTK_REAL :: dphi1,dphi2,dphi3,kdiag_bg
  CCTK_REAL :: phi_ijk,deltaijk,delta_velijk(3)
  CCTK_REAL :: delta_test,dvel_test
  CCTK_REAL, parameter :: perturb_size_tol = 1.e-5 ! perturbations larger than this tol will warn user
  CCTK_INT  :: ncells(3)
  character(len=400) :: warn_message

  call CCTK_INFO("Initialising a linearly-perturbed FLRW spacetime with a SINGLE-MODE perturbation")

  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call FLRW_SetLogicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set logicals that are specific to this single mode case only
  perturb_x = .False.; perturb_y = .False.
  perturb_z = .False.; perturb_all = .False.

  perturb_x   = CCTK_EQUALS (FLRW_perturb_direction, "x")
  perturb_y   = CCTK_EQUALS (FLRW_perturb_direction, "y")
  perturb_z   = CCTK_EQUALS (FLRW_perturb_direction, "z")
  perturb_all = CCTK_EQUALS (FLRW_perturb_direction, "all")
  !
  if (perturb_x) then
     call CCTK_INFO("    Perturbing in the x-direction ONLY ... ")
  elseif (perturb_y) then
     call CCTK_INFO("    Perturbing in the y-direction ONLY ... ")
  elseif (perturb_z) then
     call CCTK_INFO("    Perturbing in the z-direction ONLY ... ")
  elseif (perturb_all) then
     call CCTK_INFO("    Perturbing in all directions ... ")
  else
     call CCTK_INFO( "    Perturbing in no directions? " )
  endif

  !
  ! set some parameters used in setting the background
  call FLRW_SetBackground(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)

  !
  ! wavenumber is the same in each direction
  lambda = single_perturb_wavelength * boxlen
  kval   = 2.0d0 * pi / lambda(1)
  kx = kval; ky = kval; kz = kval

  ! factors for the density and velocity perturbations, respectively: eqns. (28),(29) in Macpherson+(2017)
  !    note \delta \propto kval^2 \phi not |k| because this is not technically a single mode in k-space
  !perturb_rho0 = - kval**2 / (4.0d0 * pi * rho0 * asq) - 2.0d0
  !perturb_v0   = - hub / ( 4.0d0 * pi * rhostar )
  perturb_rho0 = - 2.0d0 * kval**2 / (3.0d0 * hub**2)  - 2.0d0
  perturb_v0   = - 2.0d0 / (3.0d0 * a0 * hub)

  !
  ! rough check of the amplitude of the perturbs, and warn if larger
  !
  delta_test = abs( perturb_rho0 * phi_amplitude )
  dvel_test  = abs( perturb_v0 * kx * phi_amplitude )
  if (delta_test > perturb_size_tol .or. dvel_test > perturb_size_tol .or. phi_amplitude > perturb_size_tol) then
     write(warn_message,'(a,e12.5,a,e12.5,a,e12.5)')"Perturbation/s may be too large to satisfy linear approximation. Density:",&
          delta_test," velocity:",dvel_test," and phi:",phi_amplitude
     call CCTK_WARN(CCTK_WARN_ALERT,warn_message)
     call CCTK_WARN(CCTK_WARN_ALERT,"Please REDUCE phi_amplitude or INCREASE FLRW_init_HL to reduce perturbation sizes")
  endif

  !
  ! spatial loop over *local* grid size for this processor
  !
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !
           ! set up perturbations
           !
           phi_ijk = 0.0d0; delta_velijk = 0.0d0; deltaijk = 0.0d0

           if (perturb_x) then
              phi_ijk = phi_amplitude * sin(kx * x(i,j,k) - phi_phase_offset)
              dphi1   = phi_amplitude * kx * cos(kx * x(i,j,k) - phi_phase_offset)
              delta_velijk(1) = perturb_v0 * dphi1
           elseif (perturb_y) then
              phi_ijk = phi_amplitude * sin(ky * y(i,j,k) - phi_phase_offset)
              dphi2   = phi_amplitude * ky * cos(ky * y(i,j,k) - phi_phase_offset)
              delta_velijk(2) = perturb_v0 * dphi2
           elseif (perturb_z) then
              phi_ijk = phi_amplitude * sin(kz * z(i,j,k) - phi_phase_offset)
              dphi3   = phi_amplitude * kz * cos(kz * z(i,j,k) - phi_phase_offset)
              delta_velijk(3) = perturb_v0 * dphi3
           elseif (perturb_all) then
              phi_ijk = phi_amplitude * (sin(kx * x(i,j,k) - phi_phase_offset) + &
                   sin(ky * y(i,j,k) - phi_phase_offset) + sin(kz * z(i,j,k) - phi_phase_offset))
              dphi1   = phi_amplitude * kx * cos(kx * x(i,j,k) - phi_phase_offset)
              dphi2   = phi_amplitude * ky * cos(ky * y(i,j,k) - phi_phase_offset)
              dphi3   = phi_amplitude * kz * cos(kz * z(i,j,k) - phi_phase_offset)
              delta_velijk(1) = perturb_v0 * dphi1
              delta_velijk(2) = perturb_v0 * dphi2
              delta_velijk(3) = perturb_v0 * dphi3
           endif
           deltaijk = perturb_rho0 * phi_ijk

           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              if (lapse) then
                 alp(i,j,k) = FLRW_lapse_value * sqrt(1.0d0 + 2.0d0 * phi_ijk)
              endif

              ! time deriv of lapse -- evolution of this is specified in ADMBase.
              if (dtlapse) then
                 dtalp(i,j,k) = 0.0d0
              endif

              ! shift vector, always zero in this thorn
              if (shift) then
                 betax(i,j,k) = 0.0d0
                 betay(i,j,k) = 0.0d0
                 betaz(i,j,k) = 0.0d0
              endif

              ! set perturbed metric and K_ij
              gxx(i,j,k) = asq * (1.0d0 - 2.0d0 * phi_ijk)
              gxy(i,j,k) = 0.0d0
              gxz(i,j,k) = 0.0d0
              gyy(i,j,k) = asq * (1.0d0 - 2.0d0 * phi_ijk)
              gyz(i,j,k) = 0.0d0
              gzz(i,j,k) = asq * (1.0d0 - 2.0d0 * phi_ijk)

              kdiag_bg   = - adot * a0 / alp(i,j,k)
              kxx(i,j,k) = kdiag_bg * (1.0d0 - 2.0d0 * phi_ijk)
              kxy(i,j,k) = 0.0d0
              kxz(i,j,k) = 0.0d0
              kyy(i,j,k) = kdiag_bg * (1.0d0 - 2.0d0 * phi_ijk)
              kyz(i,j,k) = 0.0d0
              kzz(i,j,k) = kdiag_bg * (1.0d0 - 2.0d0 * phi_ijk)

              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! perturb the matter
                 press(i,j,k) = 0.0d0 ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0.0d0
                 rho(i,j,k)   = rho0 * (1.0d0 + deltaijk)
                 vel(i,j,k,:) = delta_velijk(:)
              endif

           endif

        enddo
     enddo
  enddo


end subroutine FLRW_SingleMode
