! file    init_tools.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    A collection of subroutines to do things that are used by some or all IC's in this thorn

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module FLRW_InitTools
  !use, intrinsic :: iso_c_binding    
  implicit none
  CCTK_REAL, parameter :: pi = 4.0d0 * atan(1.0d0)

contains

  subroutine FLRW_SetLogicals(lapse,dtlapse,shift,data,hydro)
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

  end subroutine FLRW_SetLogicals



  subroutine FLRW_SetBackground(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)
    !
    ! a subroutine to set common background paramters used by all (or >1) cases
    !
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL, intent(inout) :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
    CCTK_INT, intent(out) :: ncells(3)

    ncells    = cctk_gsh - 2 * cctk_nghostzones      ! Number of grid cells, resolution
    boxlen    = ncells * cctk_delta_space            ! Length of box in code units

    a0      = FLRW_init_a                         ! Initial scale factor
    asq     = a0*a0                               ! Scale factor squared
    hub     = FLRW_init_HL / boxlen(1)            ! Initial Hubble parameter
    rho0    = 3.0d0 * hub**2 / (8.0d0 * pi * asq) ! Initial background density from Friedmann eqn.
    rhostar = rho0 * a0**3                        ! Conserved FLRW density
    adot    = hub * a0                            ! a' (conformal) from H = a'/a
    hubdot  = -4.0d0 * pi * rho0 * asq / 3.0d0    ! H' from derivative of Friedmann eqns

  end subroutine FLRW_SetBackground


  !
  ! A subroutine to take in an array and do linear interpolation to a point
  !
  subroutine FLRW_Interp1DLinear(xi,nx,xvals,func,func_interp)
      integer, intent(in) :: nx        ! the size of the array
      CCTK_REAL, intent(in) :: xi ! the point of interpolation
      CCTK_REAL, intent(in) :: xvals(nx),func(nx) ! xvals must be increasing
      CCTK_REAL, intent(out) :: func_interp ! interpolated point
      CCTK_REAL :: xl,xu,xdu,xdl
      integer :: i,il,iu

      ! We assume that xi lies somewhere in the range xvals(1:)
      !    but check anyway
      if (xi<xvals(1) .or. xi>xvals(nx)) then
          call CCTK_WARN(CCTK_WARN_ALERT,"k value out of interpolation range")
      endif

      ! Find the xl and xu values to use to interpolate
      i = 1
      do while(xvals(i)<=xi)
          ! loop thru xvals until we have passed xi
          i = i + 1
      enddo
      ! note we will be 1 step *past* where we want to be, so need to -1 off indices we want
      il = i-1; iu = i
      xl = xvals(il)  ! lower limit
      xu = xvals(iu)  ! upper limit
      ! do the interpolation
      !call interp1Dlin(xi,xl,xu,(/func(il),func(iu)/),func_interp)

      !
      ! distance of interp point from edges
      xdl = (xi - xl) / (xu - xl)
      xdu = (xu - xi) / (xu - xl)
      
      ! do the interpolation
      func_interp = func(il) * xdu + func(iu) * xdl

    end subroutine FLRW_Interp1DLinear


    !
    ! A subroutine to return an array of 3D random numbers drawn from normal distribution
    !
    subroutine FLRW_GetRandomNormal3D(nx,ny,nz,randnums,rseed)
        integer, intent(in) :: nx,ny,nz,rseed
        CCTK_REAL, intent(out) :: randnums(nx,ny,nz)
        CCTK_REAL :: rand1(nx,ny,nz),rand2(nx,ny,nz)
        integer :: n
        integer, allocatable :: seedn(:)
        !
        ! 0. Initialise the random seed
        !
        ! Seed the intrinsic RANDOM_NUMBER generator using intrinsic RANDOM_SEED
        call RANDOM_SEED(size=n)
        allocate(seedn(n))
        ! apply the chosen random seed
        seedn = rseed
        call RANDOM_SEED(put=seedn)
        deallocate(seedn)
        
        !
        ! 1. Generate the random numbers in [0,1]
        !       --> when seeded, two subsequent calls give different numbers, but are the same if re-run with same seed
        call RANDOM_NUMBER(rand1)
        call RANDOM_NUMBER(rand2)
        
        !
        ! 2. Now get the normally distributed numbers using the Box-Muller transformation from a distribution in [0,1]
        randnums = sqrt(-2.0d0*log(rand1))*cos(2.0d0*pi*rand2)

      end subroutine FLRW_GetRandomNormal3D
    

end module FLRW_InitTools
