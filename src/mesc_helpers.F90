! file    mesc_helpers.F90
! author  Hayley Macpherson
! date    05.06.2021
! desc    Some routines COPIED from mescaline that we need to generate exact ICs in FLRWSolver

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module mesc_helpers
  implicit none

contains
    ! --------------------------------------------------------------------------------
    !              ROUTINES COPIED FROM MESCALINE IN JUNE 2021
    !     All routines copied have been used extensively for years within mesc,
    !       we don't expect any changes to have to be made.
    !
    ! ** CCTK updates: only variable TYPES are changed, all else is consistent w/ mesc
    ! --------------------------------------------------------------------------------


    !
    ! subroutine to calculate a given component of the (3,3) Ricci tensor
    !
    ! ** need to change to CCTK_REAL here
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
      term1 = 0._dp; term2 = 0._dp
      term3 = 0._dp; term4 = 0._dp
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
      CCTK_REAL intent(in), dimension(6,nx,nx,nx) :: gij
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

      chr = 0._dp
      do l=1,3
          call get_deriv1fourth(g3dum_p1(l),g3dum_p2(l),g3dum_m1(l),g3dum_m2(l),dx,trm1)
          call get_deriv1fourth(gdum2_p1(l),gdum2_p2(l),gdum2_m1(l),gdum2_m2(l),dx,trm2)
          call get_deriv1fourth(g23p1(l),g23p2(l),g23m1(l),g23m2(l),dx,trm3)
          chr = chr + 0.5_dp * upg1dum(l) * (trm1 + trm2 - trm3)  ! summation
      enddo

    end subroutine get_christoffel


    !
    ! function to return 4th order approx of first derivative
    !
    !real(c_double) function deriv1fourth(fp1,fp2,fm1,fm2,h)
    ! again, convert to a subroutine...
    subroutine get_deriv1fourth(fp1,fp2,fm1,fm2,h,df)
      CCTK_REAL, intent(in) :: fp1,fm1,fp2,fm2,h
      CCTK_REAL, intent(out) :: df

      df = (-fp2 + 8._dp * fp1 - 8._dp * fm1 + fm2) / (12._dp * h)

    !end function deriv1fourth
    end subroutine get_deriv1fourth


    !
    ! Calculate the *whole* (3,3) tensor A^i_j from A_{ij}
    !
    subroutine calc_up_down(Add,upgij,Aud)
      CCTK_REAL, intent(in), dimension(3,3) :: Add, upgij
      CCTK_REAL, intent(out) :: Aud(3,3)
      integer :: l,m,k
      Aud = 0._dp
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
    !real(c_double) function trace(A,gdown)
    ! converted this to a subroutine because I'm not sure about the CCTK function stuff...
    subroutine get_trace(A,gdown,traceA)
      CCTK_REAL, intent(in) :: A(6)    ! (3,3) symmetric tensor components
      CCTK_REAL, intent(in) :: gdown(3,3) !,gij(6,nx,nx,nx)
      CCTK_REAL, intent(out) :: traceA
      CCTK_REAL :: gup(3,3), detg

      !call get_metric_at_pos(i,j,k,nx,gij,gdown)
      call inv3x3(gdown,gup,detg)

      ! A = (Axx, Axy, Axz, Ayy, Ayz, Azz)
      traceA = gup(1,1) * A(1) + 2._dp * gup(1,2) * A(2) + 2._dp * gup(1,3) * A(3) + &
           & gup(2,2) * A(4) + 2._dp * gup(2,3) * A(5) + gup(3,3) * A(6)
    !end function trace
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
    ! PERIODIC BC ROUTINES COPIED FROM mescaline/src/periodic.f90
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


end module mesc_helpers
