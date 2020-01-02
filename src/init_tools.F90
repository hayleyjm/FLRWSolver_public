! file    init_tools.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    A collection of subroutines to do things that are used by all IC's in this thorn

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module init_tools

  implicit none
  integer,   parameter :: dp = 8
  CCTK_REAL, parameter :: pi = 4.*atan(1.)

contains


  subroutine set_logicals(lapse,dtlapse,shift,data,hydro)
    !
    ! a subroutine to set the logical parameters that are set by choices in the .par
    !   file that we use to create initial data in all (or >1) cases
    !
    implicit none
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
    ! we don't need info about grid, or grid variables here, so sbouldn't need to declate ARGUMENTS
    
    logical, intent(out) :: lapse,dtlapse,shift,data,hydro

    !
    ! initialise
    lapse = .False.; dtlapse = .False.; shift = .False.
    data  = .False.; hydro   = .False.

    !
    ! check what user has set in the parameter file
    lapse   = CCTK_EQUALS (initial_lapse, "flrw")
    dtlapse = CCTK_EQUALS (initial_dtlapse, "flrw")
    shift   = CCTK_EQUALS (initial_shift, "flrw")
    data    = CCTK_EQUALS (initial_data,  "flrw")
    hydro   = CCTK_EQUALS (initial_hydro, "flrw")

  end subroutine set_logicals



  subroutine set_parameters(a0,rho0,asq,rhostar,hub,adot,hubdot)
    !
    ! a subroutine to set common paramters used by all (or >1) of the initial data cases
    !
    implicit none
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    CCTK_REAL, intent(inout) :: a0,rho0,asq,rhostar,hub,adot,hubdot

    a0      = 1._dp                             ! Initial scale factor
    rho0    = FLRW_init_rho                     ! Initial background density (code units)
    asq     = a0*a0                             ! Scale factor squared
    rhostar = rho0 * a0**3                      ! Conserved FLRW density
    hub     = sqrt((8._dp * pi)*rho0*asq/3._dp) ! H (conformal) from Friedmann eqns
    adot    = hub * a0                          ! a' (conformal) from H = a'/a
    hubdot  = -4._dp * pi * rho0 * asq / 3._dp  ! H' from derivative of Friedmann eqns

  end subroutine set_parameters



  subroutine check_metric()
    !
    ! a subroutine to check the metric, something that we do at the end of each initial data routine
    !
    implicit none
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    
    if (CCTK_EQUALS (metric_type, "physical")) then
       ! do nothing
    else
       call CCTK_WARN (0, "Unknown value of ADMBase::metric_type -- FLRW only set-up for metric_type = physical")
    endif
  end subroutine check_metric
  

end module init_tools
