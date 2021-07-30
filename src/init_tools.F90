! file    init_tools.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    A collection of subroutines to do things that are used by some or all IC's in this thorn

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module init_tools
  use, intrinsic :: iso_c_binding
  implicit none
  integer,   parameter :: dp = 8
  CCTK_REAL, parameter :: pi = 4._dp * atan(1._dp)

  interface
     subroutine call_make_ics(a_init,rhostar,box_size,resol,num_ghosts,rseed) bind(c)
       import
       real(c_double), intent(in) :: a_init,rhostar,box_size
       integer(C_INT32_T), intent(in) :: resol,num_ghosts,rseed
     end subroutine call_make_ics
  end interface

contains

  subroutine set_logicals(lapse,dtlapse,shift,data,hydro)
    !
    ! a subroutine to set the logical parameters that are set by choices in the .par
    !   file that we use to create initial data in all (or >1) cases
    !
    implicit none
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
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



  subroutine set_parameters(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)
    !
    ! a subroutine to set common background paramters used by all (or >1) cases
    !
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL, intent(inout) :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
    CCTK_INT, intent(out) :: ncells(3)

    ncells  = cctk_gsh - 2 * cctk_nghostzones      ! Number of grid cells, resolution
    boxlen  = ncells * cctk_delta_space            ! Length of box in code units

    a0      = FLRW_init_a                         ! Initial scale factor
    asq     = a0*a0                               ! Scale factor squared
    hub     = FLRW_init_HL / boxlen(1)            ! Initial Hubble parameter
    rho0    = 3._dp * hub**2 / (8._dp * pi * asq) ! Initial background density from Friedmann eqn.
    rhostar = rho0 * a0**3                        ! Conserved FLRW density
    adot    = hub * a0                            ! a' (conformal) from H = a'/a
    hubdot  = -4._dp * pi * rho0 * asq / 3._dp    ! H' from derivative of Friedmann eqns

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
