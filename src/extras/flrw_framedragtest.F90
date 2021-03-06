! file    flrw_framedragtest.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    Single-mode perturbation to test the frame-dragging effect in longitudinal gauge. See Adamek,Barerra-Hinjosa,Macpherson,Mertens(+?)

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine FLRW_FramedragTest (CCTK_ARGUMENTS)
  USE init_tools
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer   :: i,j,k
  logical   :: lapse,dtlapse,shift,data,hydro
  
  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
  CCTK_REAL :: ky,b,L,cosky,cosky2,sinky,sinky2,gradH
  CCTK_REAL :: rho0denom,rho0num
  CCTK_INT  :: ncells(3)
  
  call CCTK_INFO("Initialising FRAMEDRAG TEST")
  
  ! 
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  call set_logicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set some parameters common to all routines relating to the background
  call set_parameters(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)

  !
  ! keep these to simplify things and make eqns. look as much like paper as possible
  b  = phi_amplitude
  L  = boxlen(2)
  ky = 2._dp * pi / L
  
  !
  ! spatial loop over *local* grid size for this processor
  !
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !
           ! store these cos we use them a bit in rho,vel etc
           cosky  = cos(ky * y(i,j,k))
           cosky2 = cosky * cosky
           sinky  = sin(ky * y(i,j,k))
           sinky2 = sinky * sinky
           ! set grad vector \Delta_i H_j -- this is the component of spatial metric \gamma_xy = \gamma_yx
           gradH  = (b / (hub * L)) * cosky
           
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

              ! set perturbed metric and K_ij - only grad vector in gamma_ij
              gxx(i,j,k) = asq
              gxy(i,j,k) = gradH
              gxz(i,j,k) = 0._dp
              gyy(i,j,k) = asq + gradH**2
              gyz(i,j,k) = 0._dp
              gzz(i,j,k) = asq

              kxx(i,j,k) = -adot
              kxy(i,j,k) = -b * cosky / (4._dp * L)
              kxz(i,j,k) = 0._dp
              !kyy(i,j,k) = -adot + b**2 * cosky2 / (2._dp * hub * L**2) ! Kyy from draft (hard-wired to get rho,v below)
              kyy(i,j,k) = -adot - b**2 * cosky2 / (2._dp * hub * L**2)  ! Kyy directly from RGTC
              kyz(i,j,k) = 0._dp
              kzz(i,j,k) = -adot

              !
              ! set up  matter variables
              !
              if (hydro) then

                 press(i,j,k) = 0._dp ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0._dp

                 ! set density, velocity from framedrag test (Mathematica/draft paper)
                 !
                 ! denominator and numerator in first term in notes
                 ! rho0 directly from eq.8 in draft (3.Oct.2018 -- checked 9.Mar.2020)
                 !rho0denom  = 128._dp * pi * L**2 * (16._dp * hub**2 * L**2 - 3._dp * b**2 * cosky2)
                 !rho0num    = (16._dp * hub**2 * L**2 - 3._dp * b**2 * cosky2)**2 - 64._dp * pi**2 * b**2 * sinky2
                 !rho(i,j,k) = 3._dp * rho0num / rho0denom

                 ! rho0 we get from Kyy directly from RGTC above
                 rho0num    = - 0.25_dp * ( (23._dp / 2._dp) * b**2 * cosky2 + 24._dp * L**2 * hub**2 )**2 +&
                      36._dp * b**2 * pi**2 * hub**4 * sinky2
                 rho0denom  = 16._dp * L**2 * hub**4 * pi * ( -(23._dp / 2._dp) * b**2 * cosky2 - 24._dp * L**2 * hub**2 )
                 rho(i,j,k) = rho0num / rho0denom

                 !
                 ! vel is v^i = u^i / (alp u^t) from mathematica with same u^x as draft (3.Oct.2018 -- checked 9.Mar.2020)
                 !vel(i,j,k,1) = 8._dp * pi * b * sinky / (3._dp * b**2 * cosky2 - 16._dp * L**2 * hub**2)

                 ! v^x we get from Kyy directly from RGTC above
                 vel(i,j,k,1) = -6._dp * pi * b * sinky / ((23._dp/4._dp) * b**2 * cosky2 + 12._dp * L**2 * hub**2)
                 vel(i,j,k,2) = 0._dp
                 vel(i,j,k,3) = 0._dp
              endif

           endif

        enddo
     enddo
  enddo

  !
  ! make sure the metric type is physical (not conformal)
  !
  call check_metric()


end subroutine FLRW_FramedragTest
