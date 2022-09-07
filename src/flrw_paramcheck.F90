! file    flrw_paramcheck.F90
! author  Hayley Macpherson
! date    06.09.2022
! desc    Check some parameters for thorn FLRWSolver

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine FLRW_ParamCheck (CCTK_ARGUMENTS)
  !
  ! a subroutine to check the metric type is set to physical
  !
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  if (CCTK_EQUALS (metric_type, "physical")) then
     ! do nothing
  else
     call CCTK_WARN (0, "Unknown value of ADMBase::metric_type -- FLRWSolver only set-up for metric_type = physical")
  endif
  
end subroutine FLRW_ParamCheck
