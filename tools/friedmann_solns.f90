!! program to calculate exact solutions to Friedmann equations, scale factor and density over time and 's'.
!! REQUIRES:: nothing.
program friedmann_solns
  implicit none

  integer,parameter::n=int(1e5),dp=8
  real(dp) :: a, rho, t, s, dt, res
  real(dp), parameter :: rho0=1.e-8, a0=1.,pi=3.141592653589793238462643383279502884197
  integer :: i

  print*,'Enter resolution of run (dn):'
  read(*,*) res

  dt = 0.1_dp * res
  
  open(unit=10,file='friedmann_exact_solutions.out',status='replace')
  write(10,*) 'tf        s         a         rho'
  do i=1,n
     t = t + dt*(i-1)
     s = 1._dp + sqrt(6._dp * pi * rho0) * t
     a = a0 * s**(2._dp/3._dp)
     rho = rho0 * a0**3 / (a**3)

     print*,'time is:',t
     write(10,*) t, s, a, rho

  enddo
  close(10)

end program friedmann_solns
