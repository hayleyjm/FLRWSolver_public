! file    exact_ics.F90
! author  Hayley Macpherson
! date    05.06.2021
! desc    A collection of subroutines based on the Mescaline initial conditions generator, and to call the relevant mescaline routines for exact initial data

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module exact_ics
  implicit none

contains

    !
    ! A subroutine to be called from FLRW_powerspectrum to generate exact ICs from an
    !    initially linear-generated \phi, assuming scalar-only perturbations
    !
    !  -- this subroutine has been COPIED and modified for FLRWSolver from mescaline/src/ics/mescaline_ics.f90
    !
    subroutine mescaline_ics(a0,asq,adot,phi,rhoR,vel1,vel2,vel3)
        use init_tools, only:dp
        use mesc_helpers, only:get_metric_at_pos,trace
        use contraint_violation_ICs, only:calc_mom_source,calc_ham_source
        implicit none
        DECLARE_CCTK_ARGUMENTS
        DECLARE_CCTK_FUNCTIONS
        DECLARE_CCTK_PARAMETERS

        CCTK_REAL, intent(in) :: phi(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)),&
            & a0,asq,adot

        CCTK_REAL, intent(out), dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: rhoR,vel1,vel2,vel3

        CCTK_REAL, allocatable, dimension(:,:,:,:) :: gij,kij
        CCTK_REAL, allocatable, dimension(:,:,:,:,:,:) :: Chrsijk
        CCTK_REAL, allocatable, dimension(:,:,:) :: alp,trK

        CCTK_REAL :: trRijk,rhoijk,vsq,vels(3),dx
        CCTK_REAL :: rhoAij,Add(3,3),AUU(3,3)
        CCTK_REAL :: rxx,rxy,rxz,ryy,ryz,rzz,rdd(6)
        CCTK_REAL, dimension(3,3) :: gup,gdown,kud
        CCTK_REAL :: Sd(3),Su(3),SuSd
        CCTK_REAL :: partialud(3,3),theta !! dummys - not really in use!
        CCTK_REAL :: gdown_atstencil(13,3,3),kdown_atstencil(13,3,3)

        integer :: i,j,k,l,m,n,nx,ny,nz

        !
        ! Make some things about the grid easier to access
        nx = cctk_gsh(1); ny = cctk_gsh(2); nz = cctk_gsh(3)
        dx = cctk_delta_space

        !
        ! Allocate memory
        allocate(gij(6,nx,ny,nz),kij(6,nx,ny,nz),alp(nx,ny,nz),trK(nx,ny,nz),&
            & Chrsijk(3,3,3,nx,ny,nz))

        !
        ! First, set up the metric, kij, alp, christoffels using the input phi
        !    and assuming scalar-only perturbations
        ! initialise gij,kij components
        gij = 0._dp; kij = 0._dp
        do k=1,nz
           do j=1,ny
              do i=1,nx

                 !
                 ! set up metric and curvature
                 !
                 alp(i,j,k)   = a0 * sqrt(1._dp + 2._dp * phi(i,j,k))
                 gij(1,i,j,k) = asq * (1._dp - 2._dp * phi(i,j,k))
                 gij(4,i,j,k) = asq * (1._dp - 2._dp * phi(i,j,k))
                 gij(6,i,j,k) = asq * (1._dp - 2._dp * phi(i,j,k))
                 !
                 kij(1,i,j,k) = - adot * a0 * (1._dp - 2._dp * phi(i,j,k)) / alp(i,j,k)
                 kij(4,i,j,k) = - adot * a0 * (1._dp - 2._dp * phi(i,j,k)) / alp(i,j,k)
                 kij(6,i,j,k) = - adot * a0 * (1._dp - 2._dp * phi(i,j,k)) / alp(i,j,k)

                 ! note gij here is whole grid -- but should be ok because it's just indexed
                 !     at i,j,k which we've already set
                 call get_metric_at_pos(i,j,k,nx,gij,gdown)
                 trK(i,j,k)   = trace(kij(:,i,j,k),gdown)
              enddo
           enddo
        enddo

        !
        ! Now get Christoffel symbols
        !    note we need gij at the whole grid to take derivatives
        !   and we then need Christoffels at the whole grid to get R_ij components..
        do k=1,nz
           do j=1,ny
              do i=1,nx

                 do l=1,3
                    do m=1,3
                       do n=1,3
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
        rhoR = 0._dp; vel1 = 0._dp; vel2 = 0._dp; vel3 = 0._dp
        do k=1,nz
           do j=1,ny
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

                 call get_metric_at_pos(i,j,k,ngrid,gij,gdown)
                 trRijk = trace(rdd,gdown)
                 !
                 ! get Momentum source term - Sd should be zero for synchronous
                 call calc_mom_source(i,j,k,ngrid,dx,kij,gij,gdown,Chrsijk,trK,Sd,gup,kud)
                 !
                 ! get Hamiltonian constraint violation
                 call calc_ham_source(trRijk,trK(i,j,k),kud,rhoijk)
                 !
                 ! calc S_i S^i and S^i
                 SuSd = 0._dp
                 Su   = 0._dp
                 do l=1,3
                    do m=1,3
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

        deallocate(gij,kij,alp,Chrsijk)

    end subroutine mescaline_ics





end module exact_ics
