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

  logical   :: lapse, dtlapse, shift, data, hydro
  lapse = CCTK_EQUALS (initial_lapse, "flrw")
  dtlapse = CCTK_EQUALS (initial_dtlapse, "flrw")
  shift = CCTK_EQUALS (initial_shift, "flrw")
  data  = CCTK_EQUALS (initial_data,  "flrw")
  hydro = CCTK_EQUALS (initial_hydro, "flrw")

  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
         
         !! set up metric, extrinsic curvature, lapse and shift

         if (data) then
            gxx(i,j,k) = 1.0
            gxy(i,j,k) = 0.0
            gxz(i,j,k) = 0.0
            gyy(i,j,k) = 1.0
            gyz(i,j,k) = 0.0
            gzz(i,j,k) = 1.0
            
            kxx(i,j,k) = 1.0
            kxy(i,j,k) = 0.0
            kxz(i,j,k) = 0.0
            kyy(i,j,k) = 1.0
            kyz(i,j,k) = 0.0
            kzz(i,j,k) = 1.0
         end if

         if (lapse) then
            alp(i,j,k) = 1.0
         end if
         
         ! May also need the derivative of the lapse -- this is specified in ADMBase (somehow).
         if (dtlapse) then
            dtalp(i,j,k)=1.0
         end if

         if (shift) then
            betax(i,j,k) = 0.0
            betay(i,j,k) = 0.0
            betaz(i,j,k) = 0.0
         end if
         
         ! set up  matter variables
         if (hydro) then
            rho(i,j,k) = FLRW_init_rho
            press(i,j,k) = 0.0
            eps(i,j,k) = 0.0
            vel(i,j,k,1) = 0.0
            vel(i,j,k,2) = 0.0
            vel(i,j,k,3) = 0.0
         end if
         
        end do
     end do
  end do

  if (CCTK_EQUALS (metric_type, "physical")) then
     ! do nothing
  else
     call CCTK_WARN (0, "Unknown value of ADMBase::metric_type -- FLRW only set-up for metric_type = physical")
  end if

end subroutine FLRW_InitialData
