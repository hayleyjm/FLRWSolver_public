! file    contraint_violation_ICs.F90
! author  Hayley Macpherson
! date    05.06.2021
! desc    Module to compute the matter source terms in M_i and H given a metric

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!----------------------------------------
! ** COPIED from mescaline/src/ics June 2021 - HMAC
!    edited to work with CCTK data types/ called from FLRWSolver
!----------------------------------------
module violation_ics
    use init_tools, only:dp,pi
    use mesc_helpers, only: apply_periodic,apply_periodic_fourth,&
        & get_metric_at_pos,inv3x3,get_deriv1fourth,calc_up_down
    implicit none
contains

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
       term1a = 0._dp
       ! fourth order approx
       call get_deriv1fourth(kudip1(1,l),kudip2(1,l),kudim1(1,l),kudim2(1,l),dx,dxkudxl)
       call get_deriv1fourth(kudjp1(2,l),kudjp2(2,l),kudjm1(2,l),kudjm2(2,l),dx,dykudyl)
       call get_deriv1fourth(kudkp1(3,l),kudkp2(3,l),kudkm1(3,l),kudkm2(3,l),dx,dzkudzl)
       term1a = dxkudxl + dykudyl + dzkudzl

       term1b = 0._dp; term1c = 0._dp;
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
       term2(l) = 0._dp
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
       Sd(l) = ( term1(l) - term2(l) ) / (8._dp * pi)
    enddo

  end subroutine calc_mom_source



  !
  ! subroutine to calculate Hamiltonian constraint violation IN main loop (i.e. called at every i,j,k)
  ! ** note: if sqrt() is too slow here, output both e_scale2 and Ham, calc violation in ricci.f90 at end of loop
  !
  subroutine calc_ham_source(tracer,tracek,kud,rho)
    CCTK_REAL, intent(in) :: tracer, tracek
    CCTK_REAL, intent(in), dimension(3,3) :: kud
    CCTK_REAL, intent(out) :: rho
    CCTK_REAL :: kudud
    integer :: l,m

    kudud = 0._dp
    do l=1,3
       do m=1,3
          kudud = kudud + kud(l,m) * kud(m,l)
       enddo
    enddo
    rho = ( tracer - kudud + tracek**2 ) / (16._dp*pi)

  end subroutine calc_ham_source



end module violation_ics
