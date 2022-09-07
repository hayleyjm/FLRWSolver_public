!
! file    flrw_fileread.F90
! author  Hayley Macpherson
! date    03.11.2020
! desc   Linear perturbations to FLRW initial data, longitudinal gauge, zero shift
!         read-in from supplied filenames + paths for: delta, phi, vel(3)
!
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine FLRW_FileRead (CCTK_ARGUMENTS)
  USE FLRW_InitTools
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer   :: i,j,k
  logical   :: lapse,dtlapse,shift,data,hydro
  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
  CCTK_REAL :: phi_ijk,kdiag_bg
  !
  ! globally-size arrays (to read in initial data files)
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3))   :: phi_gs, delta_gs
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3),3) :: delta_vel_gs
  !
  ! locally-sized arrays (for this processor)
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))   :: phi, delta
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),3) :: delta_vel
  !
  CCTK_INT :: ncells(3)
  character(len=200) :: deltafile,vel1file,vel2file,vel3file,phifile
  integer :: dr_unit,dv_unit1,dv_unit2,dv_unit3,p_unit
  integer :: il,jl,kl,iu,ju,ku,dlen
  
  call CCTK_INFO("Initialising linear perturbations to an FLRW spacetime")

  !
  ! Convert CCTK strings to Fortran strings for ICs filenames
  !       (send dummy length dlen since we don't use this anyway)
  !
  call CCTK_FortranString(dlen,FLRW_deltafile,deltafile)
  call CCTK_FortranString(dlen,FLRW_phifile,phifile)
  call CCTK_FortranString(dlen,FLRW_velxfile,vel1file)
  call CCTK_FortranString(dlen,FLRW_velyfile,vel2file)
  call CCTK_FortranString(dlen,FLRW_velzfile,vel3file)
  call CCTK_INFO("From files: ")
  call CCTK_INFO(trim(deltafile))
  call CCTK_INFO(trim(phifile))
  call CCTK_INFO(trim(vel1file))
  call CCTK_INFO(trim(vel2file))
  call CCTK_INFO(trim(vel3file))

  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call FLRW_SetLogicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set parameters used in setting metric, matter parameters
  ! --> note boxlen is in code units here
  !
  call FLRW_SetBackground(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)
  
  !
  ! read in perturbations from files
  !
  open(newunit=dr_unit,file=deltafile,status='old')  ! delta rho file
  open(newunit=dv_unit1,file=vel1file,status='old')  ! delta vel file [1]
  open(newunit=dv_unit2,file=vel2file,status='old')  ! delta vel file [2]
  open(newunit=dv_unit3,file=vel3file,status='old')  ! delta vel file [3] 
  open(newunit=p_unit,file=phifile,status='old')     ! phi file

  !
  ! spatial loop over *global* grid size
  do k = 1, cctk_gsh(3)
     do j = 1, cctk_gsh(2) 
        ! loop over ROW no. (i is COLUMN no.) (i,j,k) --> (column, row, z)

        !
        ! read perturbations from ICs files to global sized (gs) arrays
        !
        read(dr_unit,*)  delta_gs(:,j,k)
        read(dv_unit1,*) delta_vel_gs(:,j,k,1)
        read(dv_unit2,*) delta_vel_gs(:,j,k,2)
        read(dv_unit3,*) delta_vel_gs(:,j,k,3)
        read(p_unit,*)   phi_gs(:,j,k)
     enddo
  enddo
  call CCTK_INFO("Opened and read IC files")

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
  delta     = delta_gs(il:iu, jl:ju, kl:ku)
  delta_vel = delta_vel_gs(il:iu, jl:ju, kl:ku, :)
  phi       = phi_gs(il:iu, jl:ju, kl:ku)

  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              phi_ijk = phi(i,j,k)
              
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

              kdiag_bg   = - adot * a0 / alp(i,j,k)
              kxx(i,j,k) = kdiag_bg * (1.0d0 - 2.0d0 * phi_ijk)
              kxy(i,j,k) = 0.0d0
              kxz(i,j,k) = 0.0d0
              kyy(i,j,k) = kdiag_bg * (1.0d0 - 2.0d0 * phi_ijk)
              kyz(i,j,k) = 0.0d0
              kzz(i,j,k) = kdiag_bg * (1.0d0 - 2.0d0 * phi_ijk)

              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! perturb the matter
                 press(i,j,k) = 0.0d0 ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0.0d0
                 rho(i,j,k)   = rho0 * (1.0d0 + delta(i,j,k))
                 vel(i,j,k,:) = delta_vel(i,j,k,:)
              endif

           endif

        enddo
     enddo
  enddo


end subroutine FLRW_FileRead
