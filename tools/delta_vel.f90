program deltavel
  implicit none

  integer, parameter:: dp=8, n=int(1e5)
  integer:: i, number, j, itmin, it, l, dit, m, t, dt
  real(dp):: vel_120, vel_zero, dum, delta
  real(dp), dimension(20):: vel_x, x
  character(len=40):: file_vel

  m = 100
  dit = 512
  itmin = 0
  dt = 2.4_dp
  number = 100 !every 'number'th file

  open(unit=10, file='deltavel.out', status='replace')

  do i=1,n+1
     m = m+1
     if (i==1) then
        it = itmin
     else
        it = itmin + ((i-1)*dit*number - dit)
     endif
     t = (i-1._dp)*dt

     write(file_vel,'(a,i9.9,a)')'velocity_',it,'.dat'
     open(unit=m,file=file_vel,status='old')
     do j=1,6 ! read out headers
        read(m,*)
     enddo

     do j=1,20
        read(m,*) (dum,l=1,9), x(j), dum, dum, vel_x(j)
        if (x(j)==0) then
           vel_zero = vel_x(j)
           print*, 'vel_zero=', vel_zero
        elseif (x(j)==120) then
           vel_120 = vel_x(j)
        endif
     enddo


     delta = abs(vel_zero) - abs(vel_120)

     print*,'ITERATION is',it,'writing t=',t,'delta=',delta
     write(10,*) t, delta
     
     close(m)

  enddo

end program deltavel
     
