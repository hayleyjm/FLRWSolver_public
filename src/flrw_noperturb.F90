! file    flrw_noperturb.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    Non-perturbed FLRW initial data, longitudinal gauge, zero shift

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine FLRW_NoPerturb (CCTK_ARGUMENTS)
  USE FLRW_InitTools
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer     :: i,j,k
  logical     :: lapse,dtlapse,shift,data,hydro
  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3),kdiag_bg
  CCTK_INT  :: ncells(3)

  call CCTK_INFO("Initialising an FLRW spacetime")

  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call FLRW_SetLogicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set parameters used in setting metric, matter parameters
  !
  call FLRW_SetBackground(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)
  
  !
  ! spatial loop over *local* grid size for this processor
  !
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then
              
              if (lapse) then
                 alp(i,j,k) = FLRW_lapse_value
              endif
              
              ! time deriv of lapse -- evolution of this is specified in ADMBase.
              if (dtlapse) then
                 dtalp(i,j,k) = 0.0d0
              endif
              
              if (shift) then
                 betax(i,j,k) = 0.0d0
                 betay(i,j,k) = 0.0d0
                 betaz(i,j,k) = 0.0d0
              endif
              
              ! set FLRW metric and K_ij
              gxx(i,j,k) = asq
              gxy(i,j,k) = 0.0d0
              gxz(i,j,k) = 0.0d0
              gyy(i,j,k) = asq
              gyz(i,j,k) = 0.0d0
              gzz(i,j,k) = asq
              
              kdiag_bg   = - adot * a0 / FLRW_lapse_value
              kxx(i,j,k) = kdiag_bg
              kxy(i,j,k) = 0.0d0
              kxz(i,j,k) = 0.0d0
              kyy(i,j,k) = kdiag_bg
              kyz(i,j,k) = 0.0d0
              kzz(i,j,k) = kdiag_bg
             
              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! set matter to homogeneous values
               	 press(i,j,k) = 0.0d0 ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0.0d0
                 vel(i,j,k,:) = 0.0d0
                 rho(i,j,k)   = rho0
              endif

           endif
           
        enddo
     enddo
  enddo
  

end subroutine FLRW_NoPerturb
