!
! file    flrw_powerspectrum.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    A spectrum of perturbations to FLRW initial data, longitudinal gauge, zero shift
!            calls create_ics.py initial conditions generator for the given parameters
!
! ** a copy of the working code to test out the Fortran-only version of the code **
!
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine FLRW_Powerspectrum (CCTK_ARGUMENTS)
  USE init_tools
  USE powerspec_ics, only:make_ics
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer   :: i,j,k
  logical   :: lapse,dtlapse,shift,data,hydro
  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
  CCTK_REAL :: phi_ijk,kvalue
  !
  ! globally-size arrays (to read in initial data files)
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3))   :: phi_gs,delta_gs
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3),3) :: delta_vel_gs
  !
  ! locally-sized arrays (for this processor)
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))   :: phi, delta
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),3) :: delta_vel
  !
  CCTK_INT :: ncells(3)
  integer :: iunit,ierr
  integer :: il,jl,kl,iu,ju,ku

  call CCTK_INFO("Initialising a power spectrum of perturbations to an EdS spacetime")

  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call set_logicals(lapse,dtlapse,shift,data,hydro,callmesc)

  !
  ! set parameters used in setting metric, matter parameters
  !      --> note boxlen is in code units here
  !
  call set_parameters(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)

  !
  ! call initial conditions generator
  !      --> note FLRW_boxlength is in cMpc here
  !
  if (ncells(1)/=ncells(2)) call CCTK_WARN(CCTK_WARN_ALERT,"Non-uniform grid. We assume a uniform grid.")
  call make_ics(CCTK_ARGUMENTS,a0,hub,boxlen(1),ncells(1),delta_gs,phi_gs,delta_vel_gs(:,:,:,1),&
    & delta_vel_gs(:,:,:,2),delta_vel_gs(:,:,:,3))
  call CCTK_INFO("Done making initial conditions.")

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
  ! extract local part of global grid i.e. part this processor is using and store it
  !
  phi = phi_gs(il:iu, jl:ju, kl:ku)
  delta     = delta_gs(il:iu, jl:ju, kl:ku)
  delta_vel = delta_vel_gs(il:iu, jl:ju, kl:ku, :)
  !
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              phi_ijk = phi(i,j,k)

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
              gxx(i,j,k) = asq * (1._dp - 2._dp * phi_ijk)
              gxy(i,j,k) = 0._dp
              gxz(i,j,k) = 0._dp
              gyy(i,j,k) = asq * (1._dp - 2._dp * phi_ijk)
              gyz(i,j,k) = 0._dp
              gzz(i,j,k) = asq * (1._dp - 2._dp * phi_ijk)

              kvalue     = - adot * a0 / alp(i,j,k)
              kxx(i,j,k) = kvalue * (1._dp - 2._dp * phi_ijk)
              kxy(i,j,k) = 0._dp
              kxz(i,j,k) = 0._dp
              kyy(i,j,k) = kvalue * (1._dp - 2._dp * phi_ijk)
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
                 rho(i,j,k)   = rho0 * (1._dp + delta(i,j,k))
                 vel(i,j,k,:) = delta_vel(i,j,k,:)
              endif

           endif

        enddo
     enddo
  enddo
  !
  ! make sure the metric type is physical (not conformal)
  !
  call check_metric()


end subroutine FLRW_Powerspectrum
