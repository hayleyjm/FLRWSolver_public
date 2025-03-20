! file    flrw_exactics.F90
! author  Hayley Macpherson
! date    13.03.2025
! desc    A routine to do things we need to do to solve the constraints
!           given an initial metric perturbation \phi (called from FLRW_Powerspectrum_Exact)
!
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module FLRW_ExactICs
  implicit none

contains

    !
    ! A subroutine to be called from FLRW_powerspectrum_Exact to generate initial data given a
    !    given 3D array of phi perturbations -- then obtain \rho and v^i by solving the constraints
    !
    subroutine FLRW_SolveConstraints(CCTK_ARGUMENTS,nx,ai,asq,adot,phi,rhoR_gs,vel1_gs,vel2_gs,vel3_gs,phi_gs)
        USE FLRW_Mescaline_Helpers ! uses them all
        implicit none
        DECLARE_CCTK_ARGUMENTS
        DECLARE_CCTK_FUNCTIONS
        DECLARE_CCTK_PARAMETERS

        integer, intent(in) :: nx ! this is the size of the grid without ghosts
        CCTK_REAL, intent(in) :: phi(nx,nx,nx),ai,asq,adot

        ! global size arrays (incl ghosts) to output
        CCTK_REAL, intent(out), dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: rhoR_gs,vel1_gs,vel2_gs,vel3_gs,phi_gs
        ! global size grids without ghosts for calculations
        CCTK_REAL, dimension(nx,nx,nx) :: rhoR,vel1,vel2,vel3

        CCTK_REAL, allocatable, dimension(:,:,:,:) :: gij,kij
        CCTK_REAL, allocatable, dimension(:,:,:,:,:,:) :: Chrsijk
        CCTK_REAL, allocatable, dimension(:,:,:) :: alpha,trK
        CCTK_REAL :: trRijk,rhoijk,dx,SuSd,Su(3),Sd(3),kud(3,3)
        CCTK_REAL :: rxx,rxy,rxz,ryy,ryz,rzz,rdd(6),gdown(3,3),gup(3,3)

        integer :: i,j,k,l,m,n,ngh,istrt,ifin

        !
        ! Make some things about the grid easier to access (note nx is different to cctk_gsh)
        dx = cctk_delta_space(1)

        !
        ! Allocate memory
        allocate(gij(6,nx,nx,nx),kij(6,nx,nx,nx),alpha(nx,nx,nx),trK(nx,nx,nx),&
            & Chrsijk(3,3,3,nx,nx,nx))

        !
        ! First, set up the metric, kij, alp, christoffels using the input phi
        !    (and assuming scalar-only perturbations with dtphi=0)
        !
        gij = 0.d0; kij = 0.d0
        do k=1,nx
           do j=1,nx
              do i=1,nx

                 !
                 ! set up metric and curvature
                 !
                 alpha(i,j,k) = ai * sqrt(1.d0 + 2.d0 * phi(i,j,k))
                 gij(1,i,j,k) = asq * (1.d0 - 2.d0 * phi(i,j,k))
                 gij(4,i,j,k) = asq * (1.d0 - 2.d0 * phi(i,j,k))
                 gij(6,i,j,k) = asq * (1.d0 - 2.d0 * phi(i,j,k))
                 !
                 kij(1,i,j,k) = - adot * ai * (1.d0 - 2.d0 * phi(i,j,k)) / alpha(i,j,k)
                 kij(4,i,j,k) = - adot * ai * (1.d0 - 2.d0 * phi(i,j,k)) / alpha(i,j,k)
                 kij(6,i,j,k) = - adot * ai * (1.d0 - 2.d0 * phi(i,j,k)) / alpha(i,j,k)

                 ! note gij here is whole grid -- but should be ok because in get_metric_at_pos
                 !     it's just indexed at i,j,k which we've already set
                 call get_metric_at_pos(i,j,k,nx,gij,gdown)
                 call get_trace(kij(:,i,j,k),gdown,trK(i,j,k))
              enddo
           enddo
        enddo

        !
        ! Now get Christoffel symbols
        !    note we need gij at the whole grid to get Christoffels
        !   and we then need Christoffels at the whole grid to get R_ij components..
        !          thus this ridiculous loop is necessary
        do k=1,nx
           do j=1,nx
              do i=1,nx

                 do n=1,3
                    do m=1,3
                       do l=1,3
                          call get_christoffel(i,j,k,gij,nx,dx,l,m,n,Chrsijk(l,m,n,i,j,k))
                       enddo
                    enddo
                 enddo

              enddo
           enddo
        enddo

        !
        ! Next, get the actual matter variables from the constraint equations
        !     for the above metric, kij, gammas...
        !
        rhoR = 0.d0; vel1 = 0.d0; vel2 = 0.d0; vel3 = 0.d0
        do k=1,nx
           do j=1,nx
              do i=1,nx

                 !
                 ! calculate ricci components and calc traceR and traceK
                 call get_ricci_component(i,j,k,Chrsijk,nx,dx,1,1,rxx)
                 call get_ricci_component(i,j,k,Chrsijk,nx,dx,1,2,rxy)
                 call get_ricci_component(i,j,k,Chrsijk,nx,dx,1,3,rxz)
                 call get_ricci_component(i,j,k,Chrsijk,nx,dx,2,2,ryy)
                 call get_ricci_component(i,j,k,Chrsijk,nx,dx,2,3,ryz)
                 call get_ricci_component(i,j,k,Chrsijk,nx,dx,3,3,rzz)
                 rdd = (/ rxx, rxy, rxz, ryy, ryz, rzz /)

                 call get_metric_at_pos(i,j,k,nx,gij,gdown)
                 call get_trace(rdd,gdown,trRijk)
                 !
                 ! get Momentum density  (Sd should be zero for synchronous gauge)
                 call calc_mom_source(i,j,k,nx,dx,kij,gij,gdown,Chrsijk,trK,Sd,gup,kud)
                 !
                 ! get rho from Hamiltonian constraint
                 call calc_ham_source(trRijk,trK(i,j,k),kud,rhoijk)
                 !
                 ! calc S_i S^i and S^i
                 SuSd = 0.d0; Su   = 0.d0
                 do m=1,3
                    do l=1,3
                       SuSd = SuSd + gup(l,m) * Sd(l) * Sd(m)
                       Su(l) = Su(l) + gup(l,m) * Sd(m)
                    enddo
                 enddo
                 !
                 ! calculate rhoR from S_i terms and rho
                 rhoR(i,j,k) = rhoijk - SuSd / rhoijk
                 !
                 ! then v^i = u^i/u^t = S^i * alp / rho
                 ! BUT HydroBase v^i = u^i / alp * u^t = S^i / rho
                 vel1(i,j,k) = Su(1) / rhoijk
                 vel2(i,j,k) = Su(2) / rhoijk
                 vel3(i,j,k) = Su(3) / rhoijk
              enddo
           enddo
        enddo

        !
        ! Put data into global sized grid arrays (including boundary zones) to be passed out
        !
        ! Initialise to zero such that boundary zones remain zero
        rhoR_gs = 0.d0; vel1_gs = 0.d0; vel2_gs = 0.d0; vel3_gs = 0.d0; phi_gs = 0.d0
        ! Fill interior values
        !   Define indexing to extract interior zone
        ngh   = cctk_nghostzones(1) ! uniform grid only
        istrt = ngh+1; ifin = nx+ngh
        !  Set interior zone
        rhoR_gs(istrt:ifin,istrt:ifin,istrt:ifin) = rhoR
        vel1_gs(istrt:ifin,istrt:ifin,istrt:ifin) = vel1
        vel2_gs(istrt:ifin,istrt:ifin,istrt:ifin) = vel2
        vel3_gs(istrt:ifin,istrt:ifin,istrt:ifin) = vel3
        phi_gs(istrt:ifin,istrt:ifin,istrt:ifin)  = phi

        deallocate(gij,kij,alpha,trK,Chrsijk)

    end subroutine FLRW_SolveConstraints


end module FLRW_ExactICs
