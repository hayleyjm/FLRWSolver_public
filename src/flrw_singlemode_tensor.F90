! file    flrw_singlemode_tensor.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    Single-mode tensor ONLY perturbation of FLRW initial data, longitudinal gauge, zero shift
!            note: these ICs assume linear perturbation theory is valid --> 2nd violation of constraints

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine FLRW_SingleMode_Tensor (CCTK_ARGUMENTS)
  USE init_tools
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer   :: i,j,k
  logical   :: lapse,dtlapse,shift,data,hydro
  logical   :: perturb_x,perturb_y,perturb_z,perturb_all

  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
  CCTK_REAL :: kx,ky,kz,tnow,hsval,hsdotval
  CCTK_REAL :: hxx,hxy,hyy,dthxx,dthxy,dthyy,ksq,nfac
  CCTK_REAL :: ki(3),modk,kvalue,coskx,sinkx,lambda(3)
  CCTK_INT  :: ncells(3)
  character(len=400) :: warn_message

  call CCTK_INFO("Initialising a linearly-perturbed FLRW spacetime with a single-mode TENSOR ONLY perturbation")

  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call set_logicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set some parameters used in setting the background
  call set_parameters(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)

  !
  ! Set the wavevector based on box size. Need (k eta) < 1 for modes outside horizon
  !    -- note: this assumes the k(plus) = k(cross) (i.e. wavenumbers of 2 polarisations are equal)
  !    -- but does not assume the amplitudes of these are the same
  !
  !if (FLRW_init_HL <= 1._dp) then
    ! call CCTK_WARN(CCTK_WARN_ALERT,"Please set FLRW_init_HL > 1 for tensor perturbations outside horizon. Initial choice of hdot=0 may not be valid.")
  !endif
  lambda = boxlen              ! lambda^z only for motion in z-direction, others all zero
  kz     = 2._dp * pi / lambda(3) ! k^z for motion in only z-direction, others all zero
  ki     = 0._dp
  ki(3)  = kz
  modk   = sqrt(ki(1)**2 + ki(2)**2 + ki(3)**2)
  kx = ki(1); ky = ki(2)

  !
  ! Initial time as set by user
  tnow     = cctk_initial_time
  ! get functions hs(eta,k) and hsdot(eta,k)
  hsval    = hs(modk,tnow)
  hsdotval = hsdot(modk,tnow)

  !
  ! spatial loop over *local* grid size for this processor
  !
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

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
                 alp(i,j,k) = FLRW_lapse_value
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
              gxx(i,j,k) = asq * (1._dp + hxx)
              gxy(i,j,k) = asq * hxy
              gxz(i,j,k) = 0._dp
              gyy(i,j,k) = asq * (1._dp + hyy)
              gyz(i,j,k) = 0._dp
              gzz(i,j,k) = asq

              kvalue     = - adot * a0 / alp(i,j,k)
              kxx(i,j,k) = kvalue * (1._dp + hxx) - 0.5_dp * a0 * dthxx
              kxy(i,j,k) = kvalue * hxy - 0.5_dp * a0 * dthxy
              kxz(i,j,k) = 0._dp
              kyy(i,j,k) = kvalue * (1._dp + hyy) - 0.5_dp * a0 * dthyy
              kyz(i,j,k) = 0._dp
              kzz(i,j,k) = kvalue

              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! perturb the matter
                 press(i,j,k) = 0._dp ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0._dp
                 rho(i,j,k)   = rho0
                 vel(i,j,k,:) = 0._dp
              endif

           endif

        enddo
     enddo
  enddo

  !
  ! make sure the metric type is physical (not conformal)
  !
  call check_metric()


end subroutine FLRW_SingleMode_Tensor
