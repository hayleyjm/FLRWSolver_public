program convert_time !program to convert tlapse (out of cactus for alp=! 1) to tfriedmann to compare to exact solutions
  implicit none

  integer,parameter::n=int(1e6),dp=8
  real(dp) :: a, asq, rho, alp(n), tl(n), tf, dtl, f
  real(dp), parameter :: rho0=1.e-8, a0=1.,pi=3.141592653589793238462643383279502884197  
  integer :: i
  real(dp) :: dum, ham, mom
  tf = 0.

  open(unit=10,file='admbase::metric.maximum.asc',status='old')
  open(unit=12,file='admbase::lapse.maximum.asc',status='old')
  open(unit=13,file='hydrobase::rho.maximum.asc',status='old')
  open(unit=14,file='admconstraints::hamiltonian.norm2.asc',status='old')
  open(unit=15,file='admconstraints::momentum.norm2.asc',status='old')

  open(unit=20,file='converted_time.out',status='replace')
  write(20,*) '# lapse time,       friedmann time,          a.max,      rho.max,       lapse.max,      ham.norm2,    mom.norm2'

  do i=1,10
     read(10,*)
     read(12,*)
     read(13,*)
     read(14,*)
     read(15,*)
  enddo

  do i=1,n
     read(10,*) dum, tl(i), asq
     read(12,*) dum, dum, alp(i)
     read(13,*) dum, dum, rho
     read(14,*) dum, dum, ham
     read(15,*) dum, dum, mom
     if (i/=1) then
        dtl = tl(i) - tl(i-1)
        tf = tf + 0.5_dp * dtl * (alp(i-1) + alp(i))
     endif

     print*,'lapse time is:', tl(i), 'friedmann time is:', tf
     write(20,*) tl(i), tf, sqrt(asq), rho, alp(i), ham, mom
  enddo

end program convert_time
