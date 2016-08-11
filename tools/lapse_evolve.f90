program lapse_evolve
  implicit none

  integer, parameter:: n=int(1e5), dp=8
  real(dp),parameter :: a0=1., rho0=1.e-8, pi=3.141592653589793238462643383279502884197
  real(dp) :: a, s, alp(n), t(n), alp_exact, a3, dum, asq_num, a_num, alp_exact_num, tau, tnew, dtnew
  integer :: i,j,nt
  integer :: harmonicF, harmonicN
  tnew = 0.

  open(unit=10,file='admbase::lapse.maximum.asc',status='old')
 ! open(unit=12,file='../FLRW_alpha/admbase::metric.maximum.asc',status='old')

  print*,'Enter value of harmonicF:'
  read(*,*) harmonicF

  print*,'Enter value of harmonicN:'
  read(*,*) harmonicN

!read out headers
  do i=1,10
     read(10,*)
!     read(12,*)
  enddo

!open file for writing
  open (unit=15,file='lapse_evolution.dat',status='replace')
  write(15,*) 'time       s        tnew         lapse         lapse_exact (exp(t))'

  do i=1,n
     nt = nt + 1
     read(10,*) dum, t(i), alp(i)
!     read(12,*) dum, dum, asq_num

!     a_num = sqrt(asq_num)
     s = 1._dp + sqrt(6._dp * pi * rho0) * t(i)
!     print*,'S IS:',s

! proper time calculated from alp = dtau/dt, integrating using alp=a^3 from harmonic gauge eqn
!     tau = a0**3  * (s**3 - 1._dp) / (3._dp * sqrt(6._dp * pi * rho0))

     a = a0 * s**(2._dp/3._dp)

! calculate exact alpha based on d/dt (alp) = - harmonicF * alp**(hamonicN) * trK
     if (harmonicN==1 .and. harmonicF==-1) then
        alp_exact = log(a**(-3)) + 1._dp
     elseif (harmonicN==1 .and. harmonicF==1) then
        alp_exact = log(a**3) + 1._dp
     elseif (harmonicN==2 .and. harmonicF==-1) then
        alp_exact = a**(-3)
     elseif (harmonicN==2 .and. harmonicF==1) then
        alp_exact = a**3
     endif


     if (i/=1) then
!        do j=1,nt
        dtnew = t(i) - t(i-1)
        tnew = tnew + 0.5_dp * dtnew * (alp(i) + alp(i-1))
!        enddo
     endif

     write(15,*) t(i), s, tnew, alp(i), alp_exact
  enddo

!  tnew = 0.
!  do i=2,n
!     dtnew = t(i) - t(i-1)
!     tnew = tnew + 0.5_dp * dtnew * (alp(i) + alp(i-1))
!  enddo

end program lapse_evolve
