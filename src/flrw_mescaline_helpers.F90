! file    flrw_mescaline_helpers.F90
! author  Hayley Macpherson
! date    13.03.2025
! desc    Some routines COPIED (and AMENDED slightly) from mescaline that we need to generate exact ICs in FLRWSolver
!            [mescaline (https://github.com/hayleyjm/mescaline-1.0) is also written by H Macpherson; not plagiarism]

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module FLRW_MescalineHelpers
  use FLRW_InitTools, only: pi
  implicit none

contains

    ! --------------------------------------------------------------------------------
    !              ROUTINES COPIED FROM MESCALINE IN JUNE 2021
    !                [H Macpherson is author of mescaline]
    !          All routines are extensively tested + used since 2018
    !
    ! ** CCTK updates: only variable TYPES are changed, all else is consistent w/ mesc
    ! --------------------------------------------------------------------------------


    !
    ! subroutine to calculate a given component of the (3,3) Ricci tensor
    !
    subroutine get_ricci_component(i,j,k,Chrsijk,nx,dx,ridx1,ridx2,riccic)
      integer, intent(in) :: i,j,k,nx             ! current position in space & grid size
      CCTK_REAL, intent(in) :: Chrsijk(3,3,3,nx,nx,nx),dx
      integer, intent(in) :: ridx1, ridx2         ! the indices of the component, R_{idx1,idx2}

      CCTK_REAL, intent(out) :: riccic
      CCTK_REAL :: term1, term2, term3, term4

      integer :: l, m                             ! dummy indices of \Gamma's
      CCTK_REAL :: gamma1_p1, gamma1_m1, gamma2_p1, gamma2_m1, gamma3a, gamma3b, gamma4a, gamma4b
      CCTK_REAL :: gamma1_p2, gamma1_m2, gamma2_p2, gamma2_m2, dgamma1, dgamma2
      integer :: ip1,im1,jp1,jm1,kp1,km1,ip2,im2,jp2,jm2,kp2,km2 ! for periodic boundaries
      !
      ! A note on notation: here gamma1 indicates the christoffel symbol in the first term of the Ricci tensor (etc...):
      ! R_{idx1,idx2} = \partial_{i}\Gamma^{i}_{idx1,idx2} - \partial_{idx2}\Gamma^{l}_{idx1,l}
      ! + \Gamma^{l}_{n,l}\Gamma^{n}_{idx1,idx2} - \Gamma^{l}_{idx2,n}\Gamma^{n}_{idx1,l}
      !

      ! Apply periodic boundary conditions
      call apply_periodic(i,ip1,im1,nx)
      call apply_periodic(j,jp1,jm1,nx)
      call apply_periodic(k,kp1,km1,nx)
      call apply_periodic_fourth(i,ip2,im2,nx)
      call apply_periodic_fourth(j,jp2,jm2,nx)
      call apply_periodic_fourth(k,kp2,km2,nx)
      !
      ! initialise sums to zero
      term1 = 0.d0; term2 = 0.d0
      term3 = 0.d0; term4 = 0.d0
      do l=1,3
         if (l==1) then
            ! x-deriv: calculate \Gamma1 at i+1, i-1
            !
            ! FIRST TERM
            !call get_christoffel(ip1,j,k,gij,nx,dx,l,ridx1,ridx2,gamma1_p1)
            gamma1_p1 = Chrsijk(l,ridx1,ridx2,ip1,j,k)
            !call get_christoffel(im1,j,k,gij,nx,dx,l,ridx1,ridx2,gamma1_m1)
            gamma1_m1 = Chrsijk(l,ridx1,ridx2,im1,j,k)
            !call get_christoffel(ip2,j,k,gij,nx,dx,l,ridx1,ridx2,gamma1_p2)
            gamma1_p2 = Chrsijk(l,ridx1,ridx2,ip2,j,k)
            !call get_christoffel(im2,j,k,gij,nx,dx,l,ridx1,ridx2,gamma1_m2)
            gamma1_m2 = Chrsijk(l,ridx1,ridx2,im2,j,k)

         elseif (l==2) then
            ! y-deriv: calculate \Gamma1 at j+1, j-1
            !
            ! FIRST TERM
            !call get_christoffel(i,jp1,k,gij,nx,dx,l,ridx1,ridx2,gamma1_p1)
            gamma1_p1 = Chrsijk(l,ridx1,ridx2,i,jp1,k)
            !call get_christoffel(i,jm1,k,gij,nx,dx,l,ridx1,ridx2,gamma1_m1)
            gamma1_m1 = Chrsijk(l,ridx1,ridx2,i,jm1,k)
            !call get_christoffel(i,jp2,k,gij,nx,dx,l,ridx1,ridx2,gamma1_p2)
            gamma1_p2 = Chrsijk(l,ridx1,ridx2,i,jp2,k)
            !call get_christoffel(i,jm2,k,gij,nx,dx,l,ridx1,ridx2,gamma1_m2)
            gamma1_m2 = Chrsijk(l,ridx1,ridx2,i,jm2,k)

         elseif (l==3) then
            ! z-deriv: calculate \Gamma1 at k+1, k-1
            !
            ! FIRST TERM
            !call get_christoffel(i,j,kp1,gij,nx,dx,l,ridx1,ridx2,gamma1_p1)
            gamma1_p1 = Chrsijk(l,ridx1,ridx2,i,j,kp1)
            !call get_christoffel(i,j,km1,gij,nx,dx,l,ridx1,ridx2,gamma1_m1)
            gamma1_m1 = Chrsijk(l,ridx1,ridx2,i,j,km1)
            !call get_christoffel(i,j,kp2,gij,nx,dx,l,ridx1,ridx2,gamma1_p2)
            gamma1_p2 = Chrsijk(l,ridx1,ridx2,i,j,kp2)
            !call get_christoffel(i,j,km2,gij,nx,dx,l,ridx1,ridx2,gamma1_m2)
            gamma1_m2 = Chrsijk(l,ridx1,ridx2,i,j,km2)

         endif
         call get_deriv1fourth(gamma1_p1,gamma1_p2,gamma1_m1,gamma1_m2,dx,dgamma1)
         term1 = term1 + dgamma1

         !
         ! SECOND TERM
         if (ridx1==1) then
            !
            ! take the x-derivative of christoffel, get at i-stencil
            !call get_christoffel(ip1,j,k,gij,nx,dx,l,l,ridx2,gamma2_p1)
            gamma2_p1 = Chrsijk(l,l,ridx2,ip1,j,k)
            !call get_christoffel(im1,j,k,gij,nx,dx,l,l,ridx2,gamma2_m1)
            gamma2_m1 = Chrsijk(l,l,ridx2,im1,j,k)
            !call get_christoffel(ip2,j,k,gij,nx,dx,l,l,ridx2,gamma2_p2)
            gamma2_p2 = Chrsijk(l,l,ridx2,ip2,j,k)
            !call get_christoffel(im2,j,k,gij,nx,dx,l,l,ridx2,gamma2_m2)
            gamma2_m2 = Chrsijk(l,l,ridx2,im2,j,k)

         elseif (ridx1==2) then
            !
            ! take the y-derivaitve of christoffel, get at j-stencil
            !call get_christoffel(i,jp1,k,gij,nx,dx,l,l,ridx2,gamma2_p1)
            gamma2_p1 = Chrsijk(l,l,ridx2,i,jp1,k)
            !call get_christoffel(i,jm1,k,gij,nx,dx,l,l,ridx2,gamma2_m1)
            gamma2_m1 = Chrsijk(l,l,ridx2,i,jm1,k)
            !call get_christoffel(i,jp2,k,gij,nx,dx,l,l,ridx2,gamma2_p2)
            gamma2_p2 = Chrsijk(l,l,ridx2,i,jp2,k)
            !call get_christoffel(i,jm2,k,gij,nx,dx,l,l,ridx2,gamma2_m2)
            gamma2_m2 = Chrsijk(l,l,ridx2,i,jm2,k)

         elseif (ridx1==3) then
            !
            ! take the z-derivaitve of christoffel, get at j-stencil
            !call get_christoffel(i,j,kp1,gij,nx,dx,l,l,ridx2,gamma2_p1)
            gamma2_p1 = Chrsijk(l,l,ridx2,i,j,kp1)
            !call get_christoffel(i,j,km1,gij,nx,dx,l,l,ridx2,gamma2_m1)
            gamma2_m1 = Chrsijk(l,l,ridx2,i,j,km1)
            !call get_christoffel(i,j,kp2,gij,nx,dx,l,l,ridx2,gamma2_p2)
            gamma2_p2 = Chrsijk(l,l,ridx2,i,j,kp2)
            !call get_christoffel(i,j,km2,gij,nx,dx,l,l,ridx2,gamma2_m2)
            gamma2_m2 = Chrsijk(l,l,ridx2,i,j,km2)

         endif
         call get_deriv1fourth(gamma2_p1,gamma2_p2,gamma2_m1,gamma2_m2,dx,dgamma2)
         term2 = term2 + dgamma2

         ! additional sum for last two terms (\Gamma*\Gamma terms)
         do m=1,3
            !
            ! THIRD TERM
            !call get_christoffel(i,j,k,gij,nx,dx,l,m,l,gamma3a)
            gamma3a = Chrsijk(l,m,l,i,j,k)
            !call get_christoffel(i,j,k,gij,nx,dx,m,ridx1,ridx2,gamma3b)
            gamma3b = Chrsijk(m,ridx1,ridx2,i,j,k)

            term3 = term3 + gamma3a * gamma3b

            !
            ! FOURTH TERM
            !call get_christoffel(i,j,k,gij,nx,dx,l,ridx2,m,gamma4a)
            gamma4a = Chrsijk(l,ridx2,m,i,j,k)
            !call get_christoffel(i,j,k,gij,nx,dx,m,ridx1,l,gamma4b)
            gamma4b = Chrsijk(m,ridx1,l,i,j,k)

            term4 = term4 + gamma4a * gamma4b
         enddo
      enddo

      riccic = term1 - term2 + term3 - term4
    end subroutine get_ricci_component


    !
    ! a subroutine to take in current poisition in space, and indices of the christoffel and calc christoffel
    !
    !  --> * checked March 27th, 2020 (Hi from lockdown!) and all OK
    !
    subroutine get_christoffel(i,j,k,gij,nx,dx,idx1,idx2,idx3,chr,gupijk)
      integer, intent(in) :: nx                        ! the size of the grid (assuming uniform)
      integer, intent(in) :: i,j,k                     ! current position in space
      integer, intent(in) :: idx1, idx2, idx3          ! the indices of the christoffel: \Gamma^{idx1}_{idx2,idx3}
      CCTK_REAL, intent(in), dimension(6,nx,nx,nx) :: gij
      CCTK_REAL, intent(in) :: dx                       ! the gridspacing
      CCTK_REAL, intent(in), optional :: gupijk(3,3) ! g^ij at i,j,k -- optional because sometimes we have it, sometimes we don't
      CCTK_REAL, intent(out) :: chr                     ! the resulting christoffel symbol
      CCTK_REAL, dimension(3) :: g3dum_p1, g3dum_m1, gdum2_m1, gdum2_p1, upg1dum  ! the g_{k,l} and g_{l,j} terms we pass in depending on which derivative we are taking of these terms, see explanation of notation below...
      CCTK_REAL, dimension(3) :: g3dum_p2, g3dum_m2, gdum2_p2, gdum2_m2
      CCTK_REAL :: g23_ip1, g23_im1, g23_jp1, g23_jm1, g23_kp1, g23_km1  ! more temp gij terms
      CCTK_REAL :: g23_ip2, g23_im2, g23_jp2, g23_jm2, g23_kp2, g23_km2  ! more temp gij terms
      CCTK_REAL :: g23p1(3), g23m1(3), detg, g23p2(3), g23m2(3)
      CCTK_REAL, dimension(3,3) :: gup, gdown, gdown_ip1, gdown_im1, gdown_jp1, gdown_jm1, gdown_kp1, gdown_km1
      CCTK_REAL, dimension(3,3) :: gdown_ip2,gdown_im2,gdown_jp2,gdown_jm2,gdown_kp2,gdown_km2
      CCTK_REAL :: trm1, trm2, trm3
      integer :: l, ip1,im1,jp1,jm1,kp1,km1,ip2,im2,jp2,jm2,kp2,km2

      ! Apply periodic boundary conditions
      call apply_periodic(i,ip1,im1,nx)
      call apply_periodic(j,jp1,jm1,nx)
      call apply_periodic(k,kp1,km1,nx)
      call apply_periodic_fourth(i,ip2,im2,nx)
      call apply_periodic_fourth(j,jp2,jm2,nx)
      call apply_periodic_fourth(k,kp2,km2,nx)

      ! fill 3-metric's with values at certain positions
      !
      call get_metric_at_pos(i,j,k,nx,gij,gdown)
      call get_metric_at_pos(ip1,j,k,nx,gij,gdown_ip1)
      call get_metric_at_pos(im1,j,k,nx,gij,gdown_im1)
      call get_metric_at_pos(i,jp1,k,nx,gij,gdown_jp1)
      call get_metric_at_pos(i,jm1,k,nx,gij,gdown_jm1)
      call get_metric_at_pos(i,j,kp1,nx,gij,gdown_kp1)
      call get_metric_at_pos(i,j,km1,nx,gij,gdown_km1)
      call get_metric_at_pos(ip2,j,k,nx,gij,gdown_ip2)
      call get_metric_at_pos(im2,j,k,nx,gij,gdown_im2)
      call get_metric_at_pos(i,jp2,k,nx,gij,gdown_jp2)
      call get_metric_at_pos(i,jm2,k,nx,gij,gdown_jm2)
      call get_metric_at_pos(i,j,kp2,nx,gij,gdown_kp2)
      call get_metric_at_pos(i,j,km2,nx,gij,gdown_km2)
      call inv3x3(gdown,gup,detg)

      !
      ! Notes on notation below: g3dum is: g_{idx3,dummy}, gdum2 is: g_{dummy,idx2}
      ! where idx2, idx3 are the indices of the christoffel: \Gamma^{idx1}_{idx2, idx3}
      ! upg1dum is: g^{idx1,dummy}, g23 is: g_{idx2, idx3}
      !

      ! extract parts of the metric to send in to be calculated
      select case(idx2)
      case(1)
         ! we are taking the x-deriv of g_{k,l}
         g3dum_p1 = gdown_ip1(idx3,:)
         g3dum_m1 = gdown_im1(idx3,:)
         g3dum_p2 = gdown_ip2(idx3,:)
         g3dum_m2 = gdown_im2(idx3,:)
      case(2)
         ! we are taking the y-deriv of g_{k,l}
         g3dum_p1 = gdown_jp1(idx3,:)
         g3dum_m1 = gdown_jm1(idx3,:)
         g3dum_p2 = gdown_jp2(idx3,:)
         g3dum_m2 = gdown_jm2(idx3,:)
      case(3)
         ! we are taking the z-deriv of g_{k,l}
         g3dum_p1 = gdown_kp1(idx3,:)
         g3dum_m1 = gdown_km1(idx3,:)
         g3dum_p2 = gdown_kp2(idx3,:)
         g3dum_m2 = gdown_km2(idx3,:)
      end select

      select case(idx3)
      case(1)
         ! we are taking the x-deriv of g_{l,j}
         gdum2_p1 = gdown_ip1(:,idx2)
         gdum2_m1 = gdown_im1(:,idx2)
         gdum2_p2 = gdown_ip2(:,idx2)
         gdum2_m2 = gdown_im2(:,idx2)
       case(2)
         ! we are taking the y-deriv of g_{l,j}
         gdum2_p1 = gdown_jp1(:,idx2)
         gdum2_m1 = gdown_jm1(:,idx2)
         gdum2_p2 = gdown_jp2(:,idx2)
         gdum2_m2 = gdown_jm2(:,idx2)
      case(3)
         ! we are taking the z-deriv of g_{l,j}
         gdum2_p1 = gdown_kp1(:,idx2)
         gdum2_m1 = gdown_km1(:,idx2)
         gdum2_p2 = gdown_kp2(:,idx2)
         gdum2_m2 = gdown_km2(:,idx2)
      end select

      upg1dum = gup(idx1,:)
      g23_ip1 = gdown_ip1(idx2,idx3)
      g23_im1 = gdown_im1(idx2,idx3)
      g23_jp1 = gdown_jp1(idx2,idx3)
      g23_jm1 = gdown_jm1(idx2,idx3)
      g23_kp1 = gdown_kp1(idx2,idx3)
      g23_km1 = gdown_km1(idx2,idx3)

      g23_ip2 = gdown_ip2(idx2,idx3)
      g23_im2 = gdown_im2(idx2,idx3)
      g23_jp2 = gdown_jp2(idx2,idx3)
      g23_jm2 = gdown_jm2(idx2,idx3)
      g23_kp2 = gdown_kp2(idx2,idx3)
      g23_km2 = gdown_km2(idx2,idx3)

      ! for the term g_{idx2,idx3} put into list to take certain derivs in loop
      g23p1 = (/ g23_ip1, g23_jp1, g23_kp1 /)
      g23m1 = (/ g23_im1, g23_jm1, g23_km1 /)
      g23p2 = (/ g23_ip2, g23_jp2, g23_kp2 /)
      g23m2 = (/ g23_im2, g23_jm2, g23_km2 /)

      chr = 0.d0
      do l=1,3
          call get_deriv1fourth(g3dum_p1(l),g3dum_p2(l),g3dum_m1(l),g3dum_m2(l),dx,trm1)
          call get_deriv1fourth(gdum2_p1(l),gdum2_p2(l),gdum2_m1(l),gdum2_m2(l),dx,trm2)
          call get_deriv1fourth(g23p1(l),g23p2(l),g23m1(l),g23m2(l),dx,trm3)
          chr = chr + 0.5d0 * upg1dum(l) * (trm1 + trm2 - trm3)  ! summation
      enddo

    end subroutine get_christoffel


    !
    ! subroutine to return 4th order approx of first derivative
    !
    subroutine get_deriv1fourth(fp1,fp2,fm1,fm2,h,df)
      CCTK_REAL, intent(in) :: fp1,fm1,fp2,fm2,h
      CCTK_REAL, intent(out) :: df

      df = (-fp2 + 8.d0 * fp1 - 8.d0 * fm1 + fm2) / (12.d0 * h)

    end subroutine get_deriv1fourth


    !
    ! Calculate the *whole* (3,3) tensor A^i_j from A_{ij}
    !
    subroutine calc_up_down(Add,upgij,Aud)
      CCTK_REAL, intent(in), dimension(3,3) :: Add, upgij
      CCTK_REAL, intent(out) :: Aud(3,3)
      integer :: l,m,k
      Aud = 0.d0
      do l=1,3
         do m=1,3
            ! sum loop
            do k=1,3
               Aud(l,m) = Aud(l,m) + upgij(l,k) * Add(k,m)
            enddo
         enddo
      enddo
    end subroutine calc_up_down



    !
    ! a subroutine to return the metric (or K_ij) in (3,3) form at entire finite diff. stencil (fourth)
    !
    subroutine get_metric_at_stencil(ipos,jpos,kpos,nx,gij,gij_atstencil)
      integer, intent(in) :: ipos, jpos, kpos, nx !! positions to get metric at + grid size
      CCTK_REAL, intent(in), dimension(6,nx,nx,nx) :: gij
      CCTK_REAL, intent(out), dimension(13,3,3) :: gij_atstencil
      integer :: im2,im1,ip1,ip2,jm2,jm1,jp1,jp2,km2,km1,kp1,kp2

      ! -----------------------------------------------------------------------------------
      !
      ! NOTE: gij_atstencil(13,3,3) first dimension is:
      !        --> (ipos,jpos,kpos),im2,im1,ip1,ip2,jm2,jm1,jp1,jp2,km2,km1,kp1,kp2
      !
      !  and those with, e.g., im2 are (im2,j,k) i.e. other indices are at current position
      !
      ! -----------------------------------------------------------------------------------

      !
      ! apply periodic boundaries -- for whole stencil here
      !
      call apply_periodic(ipos,ip1,im1,nx)
      call apply_periodic(jpos,jp1,jm1,nx)
      call apply_periodic(kpos,kp1,km1,nx)
      call apply_periodic_fourth(ipos,ip2,im2,nx)
      call apply_periodic_fourth(jpos,jp2,jm2,nx)
      call apply_periodic_fourth(kpos,kp2,km2,nx)

      !
      ! get the metric / K_ij at all stencil positions
      !
      call get_metric_at_pos(ipos,jpos,kpos,nx,gij,gij_atstencil(1,:,:))

      call get_metric_at_pos(im2,jpos,kpos,nx,gij,gij_atstencil(2,:,:))
      call get_metric_at_pos(im1,jpos,kpos,nx,gij,gij_atstencil(3,:,:))
      call get_metric_at_pos(ip1,jpos,kpos,nx,gij,gij_atstencil(4,:,:))
      call get_metric_at_pos(ip2,jpos,kpos,nx,gij,gij_atstencil(5,:,:))

      call get_metric_at_pos(ipos,jm2,kpos,nx,gij,gij_atstencil(6,:,:))
      call get_metric_at_pos(ipos,jm1,kpos,nx,gij,gij_atstencil(7,:,:))
      call get_metric_at_pos(ipos,jp1,kpos,nx,gij,gij_atstencil(8,:,:))
      call get_metric_at_pos(ipos,jp2,kpos,nx,gij,gij_atstencil(9,:,:))

      call get_metric_at_pos(ipos,jpos,km2,nx,gij,gij_atstencil(10,:,:))
      call get_metric_at_pos(ipos,jpos,km1,nx,gij,gij_atstencil(11,:,:))
      call get_metric_at_pos(ipos,jpos,kp1,nx,gij,gij_atstencil(12,:,:))
      call get_metric_at_pos(ipos,jpos,kp2,nx,gij,gij_atstencil(13,:,:))


    end subroutine get_metric_at_stencil


    !
    ! subroutine to return the metric (3,3) at a given i,j,k position in space
    !
    subroutine get_metric_at_pos(ipos,jpos,kpos,nx,gij,gij_atpos)
      integer, intent(in) :: ipos, jpos, kpos, nx !! positions to get metric at + grid size
      CCTK_REAL, intent(in), dimension(6,nx,nx,nx) :: gij
      CCTK_REAL, intent(out), dimension(3,3) :: gij_atpos

      ! gij = (gxx,gxy,gxz,gyy,gyz,gzz)
      gij_atpos(1,1) = gij(1,ipos,jpos,kpos)
      gij_atpos(1,2) = gij(2,ipos,jpos,kpos)
      gij_atpos(1,3) = gij(3,ipos,jpos,kpos)
      gij_atpos(2,1) = gij(2,ipos,jpos,kpos)
      gij_atpos(2,2) = gij(4,ipos,jpos,kpos)
      gij_atpos(2,3) = gij(5,ipos,jpos,kpos)
      gij_atpos(3,1) = gij(3,ipos,jpos,kpos)
      gij_atpos(3,2) = gij(5,ipos,jpos,kpos)
      gij_atpos(3,3) = gij(6,ipos,jpos,kpos)

    end subroutine get_metric_at_pos



    !
    ! Calculate the trace of a (3,3) tensor given the *full* metric at all positions
    !        & position in space to calculate trace
    !
    subroutine get_trace(A,gdown,traceA)
      CCTK_REAL, intent(in) :: A(6)    ! (3,3) symmetric tensor components
      CCTK_REAL, intent(in) :: gdown(3,3) !,gij(6,nx,nx,nx)
      CCTK_REAL, intent(out) :: traceA
      CCTK_REAL :: gup(3,3), detg

      !call get_metric_at_pos(i,j,k,nx,gij,gdown)
      call inv3x3(gdown,gup,detg)

      ! A = (Axx, Axy, Axz, Ayy, Ayz, Azz)
      traceA = gup(1,1) * A(1) + 2.d0 * gup(1,2) * A(2) + 2.d0 * gup(1,3) * A(3) + &
           & gup(2,2) * A(4) + 2.d0 * gup(2,3) * A(5) + gup(3,3) * A(6)

    end subroutine get_trace

    !
    ! Calculate the inverse of a (3,3) matrix and output determinant
    !
    subroutine inv3x3(A,B,det)
      CCTK_REAL, intent(in), dimension(3,3) :: A
      CCTK_REAL, intent(out), dimension(3,3) :: B ! inverse matrix
      CCTK_REAL, intent(out) :: det

      det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - &
           & A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + &
           & A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

      B(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
      B(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
      B(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
      B(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
      B(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
      B(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
      B(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
      B(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
      B(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

      B(:,:) = B(:,:)/det
    end subroutine inv3x3



    !
    ! PERIODIC BC ROUTINES FROM mescaline/src/periodic.f90
    !
    subroutine apply_periodic(j,jp1,jm1,nx)
      integer, intent(in) :: j, nx
      integer, intent(inout) :: jp1, jm1
      integer :: stp

      !
      ! Set (j+1), (j-1) depending on j value, implementing periodic BC's
      !     --> now dependent on "step" value
      !
      if (j==1) then
         jp1 = j + 1
         jm1 = nx
      elseif (j==nx) then
         jp1 = 1
         jm1 = j - 1
      else
         jp1 = j + 1
         jm1 = j - 1
      endif

    end subroutine apply_periodic


    !
    ! Periodic boundaries when using fourth order derivatives
    !
    subroutine apply_periodic_fourth(j,jp2,jm2,nx)
      integer, intent(in) :: j, nx
      integer, intent(inout) :: jp2, jm2
      integer :: nxm1
      nxm1 = nx-1

      if (j==1) then
         jp2 = j + 2
         jm2 = nx - 1
      elseif (j==2) then
         jp2 = j + 2
         jm2 = nx
      elseif (j==nxm1) then
         jp2 = 1
         jm2 = j - 2
      elseif (j==nx) then
         jp2 = 2
         jm2 = j - 2
      else
         jp2 = j + 2
         jm2 = j - 2
      endif

    end subroutine apply_periodic_fourth



    ! =======================================================
    !
    !  Functions to calculate the matter terms of the constraints
    !
    !      -- not directly taken from mescaline, but built from
    !           the codes to calculate H and M_i, just swapped around
    !           to output \rho and S_i instead
    !
    ! =======================================================

    !
    ! Calculate the matter terms from constraints
    ! M_i = D_i K^j_i - D_j trK - 8piS_i = 0. (solve for S_i and hence v^i)
    ! (D_i is cov deriv associated with \gamma_ij)
    !   H = 3R - K_ij K^ij - K^2 + 16pirho = 0. (solve for \rho)
    !
    subroutine calc_mom_source(i,j,k,nx,dx,kij,gij,gdown,Chrsijk,tracek,Sd,gup,kud)
      integer, intent(in) :: nx,i,j,k
      CCTK_REAL, intent(in) :: dx
      CCTK_REAL, intent(in), dimension(nx,nx,nx) :: tracek
      CCTK_REAL, intent(in), dimension(6,nx,nx,nx) :: gij,kij
      CCTK_REAL, intent(in) :: Chrsijk(3,3,3,nx,nx,nx),gdown(3,3)
      CCTK_REAL, intent(out) :: Sd(3),kud(3,3),gup(3,3)
      ! rhor is rho in the REST frame. (rho in here is T_ab n^a n^b)

      CCTK_REAL :: term1a,term1b,term1c
      CCTK_REAL :: gamma1b_i,gamma1b_ii,gamma1b_iii,gamma1c_i,gamma1c_ii,gamma1c_iii
      CCTK_REAL :: trKip1,trKim1,trKjp1,trKjm1,trKkp1,trKkm1,detg
      CCTK_REAL :: trKip2,trKim2,trKjp2,trKjm2,trKkp2,trKkm2
      CCTK_REAL :: dxkudxl,dykudyl,dzkudzl
      CCTK_REAL, dimension(3)   :: term1,term2
      CCTK_REAL, dimension(3,3) :: kudip1,kudim1,kudjp1,kudjm1,kudkp1,kudkm1
      CCTK_REAL, dimension(3,3) :: kudip2,kudim2,kudjp2,kudjm2,kudkp2,kudkm2
      CCTK_REAL, dimension(3,3) :: gdownip1,gdownim1,gdownjp1,gdownjm1,gdownkp1,gdownkm1
      CCTK_REAL, dimension(3,3) :: gdownip2,gdownim2,gdownjp2,gdownjm2,gdownkp2,gdownkm2
      CCTK_REAL, dimension(3,3) :: gupip1,gupim1,gupjp1,gupjm1,gupkp1,gupkm1
      CCTK_REAL, dimension(3,3) :: gupip2,gupim2,gupkp2,gupkm2,gupjp2,gupjm2
      CCTK_REAL, dimension(3,3) :: kdown,kdownip1,kdownim1,kdownjp1,kdownjm1,kdownkp1,kdownkm1
      CCTK_REAL, dimension(3,3) :: kdownip2,kdownim2,kdownjp2,kdownjm2,kdownkp2,kdownkm2
      integer :: ip1,jp1,kp1,im1,jm1,km1,l,m
      integer :: ip2,im2,jp2,jm2,kp2,km2

      call apply_periodic(i,ip1,im1,nx)
      call apply_periodic(j,jp1,jm1,nx)
      call apply_periodic(k,kp1,km1,nx)
      call apply_periodic_fourth(i,ip2,im2,nx)
      call apply_periodic_fourth(j,jp2,jm2,nx)
      call apply_periodic_fourth(k,kp2,km2,nx)

      ! get K_ij at all positions
      call get_metric_at_pos(i,j,k,nx,kij,kdown)
      call get_metric_at_pos(ip1,j,k,nx,kij,kdownip1)
      call get_metric_at_pos(im1,j,k,nx,kij,kdownim1)
      call get_metric_at_pos(i,jp1,k,nx,kij,kdownjp1)
      call get_metric_at_pos(i,jm1,k,nx,kij,kdownjm1)
      call get_metric_at_pos(i,j,kp1,nx,kij,kdownkp1)
      call get_metric_at_pos(i,j,km1,nx,kij,kdownkm1)

      call get_metric_at_pos(ip2,j,k,nx,kij,kdownip2)
      call get_metric_at_pos(im2,j,k,nx,kij,kdownim2)
      call get_metric_at_pos(i,jp2,k,nx,kij,kdownjp2)
      call get_metric_at_pos(i,jm2,k,nx,kij,kdownjm2)
      call get_metric_at_pos(i,j,kp2,nx,kij,kdownkp2)
      call get_metric_at_pos(i,j,km2,nx,kij,kdownkm2)

      ! get g_ij and g^ij at all positions
      !call get_metric_at_pos(i,j,k,nx,gij,gdown)
      call get_metric_at_pos(ip1,j,k,nx,gij,gdownip1)
      call get_metric_at_pos(im1,j,k,nx,gij,gdownim1)
      call get_metric_at_pos(i,jp1,k,nx,gij,gdownjp1)
      call get_metric_at_pos(i,jm1,k,nx,gij,gdownjm1)
      call get_metric_at_pos(i,j,kp1,nx,gij,gdownkp1)
      call get_metric_at_pos(i,j,km1,nx,gij,gdownkm1)

      call get_metric_at_pos(ip2,j,k,nx,gij,gdownip2)
      call get_metric_at_pos(im2,j,k,nx,gij,gdownim2)
      call get_metric_at_pos(i,jp2,k,nx,gij,gdownjp2)
      call get_metric_at_pos(i,jm2,k,nx,gij,gdownjm2)
      call get_metric_at_pos(i,j,kp2,nx,gij,gdownkp2)
      call get_metric_at_pos(i,j,km2,nx,gij,gdownkm2)

      call inv3x3(gdown,gup,detg)
      call inv3x3(gdownip1,gupip1,detg)
      call inv3x3(gdownim1,gupim1,detg)
      call inv3x3(gdownjp1,gupjp1,detg)
      call inv3x3(gdownjm1,gupjm1,detg)
      call inv3x3(gdownkp1,gupkp1,detg)
      call inv3x3(gdownkm1,gupkm1,detg)

      call inv3x3(gdownip2,gupip2,detg)
      call inv3x3(gdownim2,gupim2,detg)
      call inv3x3(gdownjp2,gupjp2,detg)
      call inv3x3(gdownjm2,gupjm2,detg)
      call inv3x3(gdownkp2,gupkp2,detg)
      call inv3x3(gdownkm2,gupkm2,detg)

      ! get K^i_j at all positions
      call calc_up_down(kdown,gup,kud)
      call calc_up_down(kdownip1,gupip1,kudip1)
      call calc_up_down(kdownim1,gupim1,kudim1)
      call calc_up_down(kdownjp1,gupjp1,kudjp1)
      call calc_up_down(kdownjm1,gupjm1,kudjm1)
      call calc_up_down(kdownkp1,gupkp1,kudkp1)
      call calc_up_down(kdownkm1,gupkm1,kudkm1)

      call calc_up_down(kdownip2,gupip2,kudip2)
      call calc_up_down(kdownim2,gupim2,kudim2)
      call calc_up_down(kdownjp2,gupjp2,kudjp2)
      call calc_up_down(kdownjm2,gupjm2,kudjm2)
      call calc_up_down(kdownkp2,gupkp2,kudkp2)
      call calc_up_down(kdownkm2,gupkm2,kudkm2)

      ! get traceK at all positions
      trKip1 = tracek(ip1,j,k)
      trKim1 = tracek(im1,j,k)
      trKjp1 = tracek(i,jp1,k)
      trKjm1 = tracek(i,jm1,k)
      trKkp1 = tracek(i,j,kp1)
      trKkm1 = tracek(i,j,km1)

      trKip2 = tracek(ip2,j,k)
      trKim2 = tracek(im2,j,k)
      trKjp2 = tracek(i,jp2,k)
      trKjm2 = tracek(i,jm2,k)
      trKkp2 = tracek(i,j,kp2)
      trKkm2 = tracek(i,j,km2)

      ! loop over COMPONENTS of M_i
      do l=1,3

         ! D_j K^j_i term
         term1a = 0.d0
         ! fourth order approx
         call get_deriv1fourth(kudip1(1,l),kudip2(1,l),kudim1(1,l),kudim2(1,l),dx,dxkudxl)
         call get_deriv1fourth(kudjp1(2,l),kudjp2(2,l),kudjm1(2,l),kudjm2(2,l),dx,dykudyl)
         call get_deriv1fourth(kudkp1(3,l),kudkp2(3,l),kudkm1(3,l),kudkm2(3,l),dx,dzkudzl)
         term1a = dxkudxl + dykudyl + dzkudzl

         term1b = 0.d0; term1c = 0.d0;
         do m=1,3
            ! loop for sums in christoffel terms
            !call get_christoffel(i,j,k,gij,nx,dx,m,1,m,gamma1b_i)
            gamma1b_i = Chrsijk(m,1,m,i,j,k)
            !call get_christoffel(i,j,k,gij,nx,dx,m,2,m,gamma1b_ii)
            gamma1b_ii = Chrsijk(m,2,m,i,j,k)
            !call get_christoffel(i,j,k,gij,nx,dx,m,3,m,gamma1b_iii)
            gamma1b_iii = Chrsijk(m,3,m,i,j,k)
            term1b = term1b + ( gamma1b_i * kud(1,l) + gamma1b_ii * kud(2,l) + gamma1b_iii * kud(3,l) )

            !call get_christoffel(i,j,k,gij,nx,dx,1,l,m,gamma1c_i)
            gamma1c_i = Chrsijk(1,l,m,i,j,k)
            !call get_christoffel(i,j,k,gij,nx,dx,2,l,m,gamma1c_ii)
            gamma1c_ii = Chrsijk(2,l,m,i,j,k)
            !call get_christoffel(i,j,k,gij,nx,dx,3,l,m,gamma1c_iii)
            gamma1c_iii = Chrsijk(3,l,m,i,j,k)
            term1c = term1c + ( gamma1c_i * kud(m,1) + gamma1c_ii * kud(m,2) + gamma1c_iii * kud(m,3) )

         enddo
         term1(l) = term1a + term1b - term1c

         ! D_i trK term
         term2(l) = 0.d0
         select case(l)
         case(1)
             ! x-deriv of trK
             call get_deriv1fourth(trKip1,trKip2,trKim1,trKim2,dx,term2(l))
         case(2)
             ! y-deriv of trK
             call get_deriv1fourth(trKjp1,trKjp2,trKjm1,trKjm2,dx,term2(l))
         case(3)
             ! z-deriv of trK
             call get_deriv1fourth(trKkp1,trKkp2,trKkm1,trKkm2,dx,term2(l))
         end select
         !
         ! S_i for M_i = 0
         !
         Sd(l) = ( term1(l) - term2(l) ) / (8.d0 * pi)
      enddo

    end subroutine calc_mom_source



    !
    ! subroutine to calculate rho from Hamiltonian constraint inside a loop
    !
    subroutine calc_ham_source(tracer,tracek,kud,rho)
      CCTK_REAL, intent(in) :: tracer, tracek
      CCTK_REAL, intent(in), dimension(3,3) :: kud
      CCTK_REAL, intent(out) :: rho
      CCTK_REAL :: kudud
      integer :: l,m

      kudud = 0.d0
      do l=1,3
         do m=1,3
            kudud = kudud + kud(l,m) * kud(m,l)
         enddo
      enddo
      rho = ( tracer - kudud + tracek**2 ) / (16.d0*pi)

    end subroutine calc_ham_source


end module FLRW_MescalineHelpers
