! file    FLRW_SingleMode_Scalar_Tensor.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    Single-mode perturbation of FLRW initial data, longitudinal gauge, zero shift
!         perturbation in BOTH scalar and tensor modes
!            note: these ICs assume linear perturbation theory is valid --> O(10^-10) violation of constraints

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine FLRW_SingleMode_Scalar_Tensor (CCTK_ARGUMENTS)
  USE init_tools
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer   :: i,j,k
  logical   :: lapse,dtlapse,shift,data,hydro
  logical   :: perturb_x,perturb_y,perturb_z,perturb_all

  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
  CCTK_REAL :: lambda(3),kx,ky,kz,modk,perturb_rho0,perturb_v0
  CCTK_REAL :: nfac,hxx,hxy,hyy,tnow,hsval,hsdotval
  CCTK_REAL :: dphi1,dphi2,dphi3,kvalue,coskx,sinkx,lambda_hij(3)
  CCTK_REAL :: kx_hij,ky_hij,kz_hij,modk_hij,ksq_hij,ki_hij(3)
  CCTK_REAL :: phi_ijk,deltaijk,delta_velijk(3)
  CCTK_REAL :: delta_test,dvel_test
  CCTK_REAL, parameter :: perturb_size_tol = 1.e-5 ! perturbations larger than this tol will warn user
  CCTK_INT  :: ncells(3)
  character(len=400) :: warn_message

  call CCTK_INFO("Initialising a linearly-perturbed FLRW spacetime with a SINGLE-MODE perturbation")

  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call set_logicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set logicals that are specific to this single mode case only
  !       --> SCALAR only
  !
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
  call set_parameters(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)

  !
  ! Set the tensor wavevector based on box size
  !    -- note: this assumes the k(plus) = k(cross) (i.e. wavenumbers of 2 polarisations are equal)
  !    -- but does not assume the amplitudes of these are the same
  !
  lambda_hij = boxlen                  ! lambda^z only for motion in z-direction, others all zero
  kz_hij     = 2._dp * pi / lambda_hij(3) ! k^z for motion in only z-direction, others all zero
  ki_hij     = 0._dp
  ki_hij(3)  = kz_hij
  modk_hij   = sqrt(ki_hij(1)**2 + ki_hij(2)**2 + ki_hij(3)**2)
  kx_hij     = ki_hij(1); ky_hij = ki_hij(2)
  !
  ! Initial time as set by user
  tnow     = cctk_initial_time
  ! get functions hs(eta,k) and hsdot(eta,k) for h_ij and h'_ij
  hsval    = hs(modk_hij,tnow)
  hsdotval = hsdot(modk_hij,tnow)

  !
  ! wavenumbers for phi perturbation in each direction
  lambda = single_perturb_wavelength * boxlen
  kx = 2._dp * pi / lambda(1)
  ky = 2._dp * pi / lambda(2)
  kz = 2._dp * pi / lambda(3)
  modk = sqrt(kx**2 + ky**2 + kz**2)

  ! factors for the density and velocity perturbations, respectively: eqns. (28),(29) in Macpherson+(2017)
  perturb_rho0 = - kx**2 / (4._dp * pi * rho0 * asq) - 2._dp
  perturb_v0   = - sqrt(a0 / ( 6._dp * pi * rhostar )) / a0
  !
  ! rough check of the amplitude of scalar perturbs, and warn if larger
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
           phi_ijk = 0._dp; delta_velijk = 0._dp; deltaijk = 0._dp

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
           ! Set components of h_ij and dt_hij
           !
           coskx = cos(kx * x(i,j,k) + ky * y(i,j,k) + kz * z(i,j,k))
           sinkx = sin(kx * x(i,j,k) + ky * y(i,j,k) + kz * z(i,j,k))
           !endif
           hxx = 0.5_dp * hsval * (hplus_amplitude + hcross_amplitude) * coskx
           hyy = - hxx
           hxy = 0.5_dp * hsval * (hcross_amplitude - hplus_amplitude) * sinkx

           dthxx = 0.5_dp * hsdotval * (hplus_amplitude + hcross_amplitude) * coskx
           dthyy = - dthxx
           dthxy = 0.5_dp * hsdotval * (hcross_amplitude - hplus_amplitude) * sinkx

           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              if (lapse) then
                 alp(i,j,k) = FLRW_lapse_value * sqrt(1._dp + 2._dp * phi_ijk)
              endif

              ! time deriv of lapse -- evolution of this is specified in ADMBase.
              if (dtlapse) then
                 dtalp(i,j,k) = 0._dp
              endif

              ! shift vector, always zero in this thorn
              if (shift) then
                 betax(i,j,k) = 0._dp
                 betay(i,j,k) = 0._dp
                 betaz(i,j,k) = 0._dp
              endif

              ! set perturbed metric and K_ij
              gxx(i,j,k) = asq * (1._dp + hxx - 2._dp * phi_ijk)
              gxy(i,j,k) = asq * hxy
              gxz(i,j,k) = 0._dp
              gyy(i,j,k) = asq * (1._dp + hyy - 2._dp * phi_ijk)
              gyz(i,j,k) = 0._dp
              gzz(i,j,k) = asq * (1._dp - 2._dp * phi_ijk)

              kvalue     = - adot * a0 / alp(i,j,k)
              kxx(i,j,k) = kvalue * (1._dp + hxx - 2._dp * phi_ijk) - 0.5_dp * a0 * dthxx
              kxy(i,j,k) = kvalue * hxy - 0.5_dp * a0 * dthxy
              kxz(i,j,k) = 0._dp
              kyy(i,j,k) = kvalue * (1._dp + hyy - 2._dp * phi_ijk) - 0.5_dp * a0 * dthyy
              kyz(i,j,k) = 0._dp
              kzz(i,j,k) = kvalue * (1._dp - 2._dp * phi_ijk)

              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! perturb the matter
                 press(i,j,k) = 0._dp ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0._dp
                 rho(i,j,k)   = rho0 * (1._dp + deltaijk)
                 vel(i,j,k,:) = delta_velijk(:)
              endif

           endif

        enddo
     enddo
  enddo

  !
  ! make sure the metric type is physical (not conformal)
  !
  call check_metric()


end subroutine FLRW_SingleMode_Scalar_Tensor
