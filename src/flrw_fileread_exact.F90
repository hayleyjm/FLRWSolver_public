!
! file    flrw_fileread_exact.F90
! author  Hayley Macpherson
! date    13.03.2025
! desc    Perturbations to FLRW initial data read from input file, longitudinal gauge, zero shift
!            solves the constraints exactly given a set of metric fluctuations
!
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine FLRW_FileRead_Exact (CCTK_ARGUMENTS)
  USE FLRW_InitTools, only:FLRW_SetLogicals,FLRW_SetBackground,FLRW_Interp1DLinear
  USE FLRW_ExactICs, only:FLRW_SolveConstraints
  USE FLRW_PowerspecICs, only:FLRW_MakePkICs_Phi
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer   :: i,j,k
  logical   :: lapse,dtlapse,shift,data,hydro
  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
  CCTK_REAL :: phi_ijk
  !
  ! globally-size arrays (to read in initial data files)
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: phi_gs,rhoR_gs
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: vel1_gs,vel2_gs,vel3_gs
  !
  ! locally-sized arrays (for this processor)
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: phi,rhoR
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: vel1,vel2,vel3
  !
  ! we need a phi array that is just nx,nx,nx which is then placed into the cctk_gsh arrays
  CCTK_REAL, dimension(:,:,:), allocatable :: phi_nx
  !
  CCTK_INT :: ncells(3)
  integer :: il,jl,kl,iu,ju,ku
  integer :: dlen,ierr,p_unit
  character(len=200) :: phifile

  call CCTK_INFO("Initialising your supplied phi perturbation to an EdS spacetime")

  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call FLRW_SetLogicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set parameters used in setting metric, matter parameters
  !      --> note boxlen is in code units here
  !
  call FLRW_SetBackground(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)

  if (ncells(1)/=ncells(2)) call CCTK_WARN(CCTK_WARN_ALERT,"Non-uniform grid. We assume a uniform grid.")

  !
  ! First: set up & read in the user-supplied file specifying \phi
  call CCTK_FortranString(dlen,FLRW_phifile,phifile)
  call CCTK_INFO(trim(phifile))
  open(newunit=p_unit,file=phifile,status='old',iostat=ierr)     ! phi file
  if (ierr/=0) call CCTK_WARN(CCTK_WARN_ABORT,"Problem reading in phi perturbation")
  !
  ! spatial loop over *global* grid size
  allocate(phi_nx(ncells(1),ncells(1),ncells(1)))
  do k = 1, ncells(1)
     do j = 1, ncells(1)
        ! loop over ROW no. (i is COLUMN no.) (i,j,k) --> (column, row, z)
        read(p_unit,*)   phi_nx(:,j,k)
     enddo
  enddo
  call CCTK_INFO("Opened and read initial phi")


  !
  ! Next we need to solve the constraints to get a rho and v^i from the perturbed metric
  call FLRW_SolveConstraints(CCTK_ARGUMENTS,ncells(1),a0,asq,adot,phi_nx,rhoR_gs,vel1_gs,vel2_gs,vel3_gs,phi_gs)
  call CCTK_INFO("Done making initial conditions.")
  deallocate(phi_nx) ! don't need this no mo

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
  phi  = phi_gs(il:iu, jl:ju, kl:ku)
  rhoR = rhoR_gs(il:iu, jl:ju, kl:ku)
  vel1 = vel1_gs(il:iu, jl:ju, kl:ku)
  vel2 = vel2_gs(il:iu, jl:ju, kl:ku)
  vel3 = vel3_gs(il:iu, jl:ju, kl:ku)
  !
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              phi_ijk = phi(i,j,k)

              ! initialise the lapse given the background lapse value FLRW_lapse_value
              !      -- note we require phi_ijk << 1 so that 1+2phi>0 always
              if (lapse) then
                 alp(i,j,k) = FLRW_lapse_value * sqrt(1.0d0 + 2.0d0 * phi_ijk)
              endif

              ! time deriv of lapse -- evolution of this is specified in ADMBase.
              if (dtlapse) then
                 dtalp(i,j,k) = 0.0d0
              endif

              ! shift vector, always zero in this thorn
              if (shift) then
                 betax(i,j,k) = 0.0d0
                 betay(i,j,k) = 0.0d0
                 betaz(i,j,k) = 0.0d0
              endif

              ! set perturbed metric and K_ij
              gxx(i,j,k) = asq * (1.0d0 - 2.0d0 * phi_ijk)
              gxy(i,j,k) = 0.0d0
              gxz(i,j,k) = 0.0d0
              gyy(i,j,k) = asq * (1.0d0 - 2.0d0 * phi_ijk)
              gyz(i,j,k) = 0.0d0
              gzz(i,j,k) = asq * (1.0d0 - 2.0d0 * phi_ijk)

              kxx(i,j,k) = - adot * a0 * (1.0d0 - 2.0d0 * phi_ijk) / alp(i,j,k)
              kxy(i,j,k) = 0.0d0
              kxz(i,j,k) = 0.0d0
              kyy(i,j,k) = - adot * a0 * (1.0d0 - 2.0d0 * phi_ijk) / alp(i,j,k)
              kyz(i,j,k) = 0.0d0
              kzz(i,j,k) = - adot * a0 * (1.0d0 - 2.0d0 * phi_ijk) / alp(i,j,k)

              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! perturb the matter
                 press(i,j,k) = 0.0d0 ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0.0d0
                 rho(i,j,k)   = rhoR(i,j,k)
                 vel(i,j,k,1) = vel1(i,j,k)
                 vel(i,j,k,2) = vel2(i,j,k)
                 vel(i,j,k,3) = vel3(i,j,k)
              endif

           endif

        enddo
     enddo
  enddo


end subroutine FLRW_FileRead_Exact
