!program to convert 't_lapse' (out of cactus for dt(alp)\=0) to 't_friedmann' to compare to exact solutions. writes both times, gxx, rho, alp to output file
!! REQUIRES: --> metric.max, lapse.max, rho.max
program convert_time
  implicit none

  integer,parameter::n=int(1e5),dp=8
  real(dp) :: gxx, rho, alp(n), tl(n), tf, dtl
  real(dp), parameter :: rho0=1.e-8, a0=1.,pi=3.141592653589793238462643383279502884197  
  integer :: i
  real(dp) :: dum
  tf = 0.

  open(unit=10,file='admbase::metric.maximum.asc',status='old')
  open(unit=12,file='admbase::lapse.maximum.asc',status='old')
  open(unit=13,file='hydrobase::rho.maximum.asc',status='old')

  open(unit=20,file='converted_time.out',status='replace')
  write(20,*) 'lapse time              friedmann time            gxx              rho              alpha'

  do i=1,10
     read(10,*)
     read(12,*)
     read(13,*)
  enddo

  do i=1,n
     read(10,*) dum, tl(i), gxx
     read(12,*) dum, dum, alp(i)
     read(13,*) dum, dum, rho

     if (i/=1) then
        dtl = tl(i) - tl(i-1)
        tf = tf + 0.5_dp * dtl * (alp(i-1) + alp(i))
     endif

     print*,'lapse time is:', tl(i), 'friedmann time is:', tf
     write(20,*) tl(i), tf, gxx, rho, alp(i)
  enddo

end program convert_time
