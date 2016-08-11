! to reconstruct phi as a function of x
program phi
  implicit none

  integer, parameter :: dp=8, n=int(1e5)!41599 ! n total number of its in files
  integer :: xpts, res
  integer, parameter :: xmin=-240, xmax=240
  real(dp), allocatable, dimension(:) :: gxx_x, phi_curv, phi_gxx, x, curv_x, phi_curv_cut, phi_gxx_cut, rho_x, delta
  real(dp), allocatable, dimension(:) :: ham_x, mom_x, vel_x
  real(dp) :: dum, a, zeta, adot, t, dt, tmin, s, dads, delta_phi_curv, delta_phi_gxx, delta_rho
  real(dp) :: phi_curv_avg, phi_gxx_avg, sum, phi_curv_diff, phi_gxx_diff, delta_exact, vel_exact, curv_exact
  real(dp), parameter :: rho0=1.e-8, pi=3.141592653589793238462643383279502884197, a0=1.0, dsdt=sqrt(6.0*pi*rho0), kx=2.*pi/480.
  real(dp),parameter :: phi_exact=1.e-8
  integer :: q, p, i, j, l, m, it, itmin, dit, neg240, pos240, zero, neg120, pos120, number
  character(len=40) :: filename_metric, filename_curv, file1, file2, file3, filename_density, file_ham, file_mom, file_vel

  print*,'Enter resolution of run (dx):'
  read(*,*) res

  print*,'Every iteration(1) or every hundredth(100)?'
  read(*,*) number

  !allocate number of xpts relative to resolution
  xpts = (xmax - xmin)/res + 1
  allocate(gxx_x(xpts),phi_curv(xpts),phi_gxx(xpts),x(xpts),curv_x(xpts))
  allocate(rho_x(xpts),delta(xpts),ham_x(xpts),mom_x(xpts),vel_x(xpts))
  !arrays with boundaries cut out
  allocate(phi_curv_cut(xpts-8),phi_gxx_cut(xpts-8))

  dit = 512 !difference between iterations
  itmin = 0
  m = 100
  p = 50000
  tmin = 0.
!  dt = 2.4 ! same as cactus
  dt = 0.1 * res

  !open file for maxvals - to compare to what we see
!  open(unit=1,file='phixneg120_time.out',status='replace')
!  open(unit=2,file='phix0_time.out',status='replace')
!  open(unit=3,file='phix120_time.out',status='replace')

  ! open files to write phi(s) at each grid point
  q = 100000
  do i=1,xpts
     x(i) = xmin + (i-1)*res
!     print*,x(i)
     if (x(i)<0) then
        write(file3,'(a,i3.3,a)')'phix_neg',int(abs(x(i))),'_time.out'
 !       print*,file3
     else
        write(file3,'(a,i3.3,a)')'phix_',int(x(i)),'_time.out'
 !       print*,file3
     endif
     open(unit=q,file=file3)
     q = q+1
  enddo

  !file to look at growth of boundary values
  open (unit=100,file='phi_growth.out')
  open(unit=101,file='phi_avg.out')
  open(unit=4,file='exact_vals.out')
  open(unit=10,file='constraints_x0.out')

  do i=1,n+1
     m = m+6
     p = p+3
     it = itmin + (i-1)*dit
     t = tmin + (i-1._dp)*dt
     s = 1.0_dp + sqrt(6.0_dp * pi * rho0) * t
!     it = it + (i-1)*dit

     if (mod(i,number)==0 .or. i==2) then
 !       it = it + (i-1)*dit
        print*,'ITERATION is:',it
        if (abs(t-(1918._dp))<tiny(0._dp)) then
           print*,'iteration at t=1918=',it
        elseif (abs(t-(0._dp))<tiny(37440._dp)) then
           print*,'iteration at t=37440=',it
        elseif (abs(t-(0._dp))<tiny(11280._dp)) then
           print*,'iteration at 11280=',it
        endif

        write(filename_metric,'(a,i9.9,a)')'metric_',it,'.dat'
        open(unit=m,file=filename_metric,status='old')
        write(filename_curv,'(a,i9.9,a)')'curv_',it,'.dat'
        open(unit=m+1,file=filename_curv,status='old')
        write(filename_density,'(a,i9.9,a)')'density_',it,'.dat'
        open(unit=m+2,file=filename_density,status='old')
        write(file_ham,'(a,i9.9,a)')'hamiltonian_',it,'.dat'
        open(unit=m+3,file=file_ham,status='old')
        write(file_mom,'(a,i9.9,a)')'momentum_',it,'.dat'
        open(unit=m+4,file=file_mom,status='old')
        write(file_vel,'(a,i9.9,a)')'velocity_',it,'.dat'
        open(unit=m+5,file=file_vel,status='old')


 !    print*,'opening ',filename_metric,'and ',filename_curv, filename_density

     !read out header files
        do j=1,6
           read(m,*)
           read(m+1,*)
           read(m+2,*)
           read(m+3,*)
           read(m+4,*)
           read(m+5,*)
        enddo

        a = a0 * s**(2._dp/3._dp)
        dads = (2._dp/3._dp) * a0 * s**(-1._dp/3._dp)
        adot = dads * dsdt !chain rule for da/dt
        
        write(file1,'(a,i9.9,a)')'phicurv_',it,'.dat'
        write(file2,'(a,i9.9,a)')'phigxx_',it,'.dat'
        write(file3,'(a,i9.9,a)')'deltarho_',it,'.dat'
        open(unit=p,file=file1)
        open(unit=p+1,file=file2)
        open(unit=p+2,file=file3)
        !timestamp
        write(p,*) s
        write(p+1,*) s
        write(p+2,*) t

        q = 100000
        do j=1,xpts
           read(m,*) (dum,l=1,9), x(j), dum, dum, gxx_x(j)
           read(m+1,*) (dum,l=1,12), curv_x(j)
           read(m+2,*) (dum,l=1,12), rho_x(j)
           read(m+3,*) (dum,l=1,12), ham_x(j)
           read(m+4,*) (dum,l=1,12), mom_x(j)
           read(m+5,*) (dum,l=1,12), vel_x(j)
           !  if (i==1) print*,j,x(j),gxx_x(j)
!           print*,rho_x(j)
           zeta = (curv_x(j) / (adot * a))**2
           phi_curv(j) = 0.25_dp * (2.0_dp + zeta - sqrt(8.0_dp * zeta + zeta**2))
           phi_gxx(j) = 0.5_dp * (1.0_dp - (gxx_x(j) / (a**2) ) )
           
           do l=1,xpts
              if (l==j) then
                 write(q,*) s, x(j), phi_curv(j)
                 !              print*, s, x(j), phi_curv(j)
                 q = q+1
              endif
           enddo
           
           write(p,*) x(j), phi_curv(j)
           write(p+1,*) x(j), phi_gxx(j)
        enddo
        !end of x-loop
        
        !cut out boundaries
        phi_curv_cut = phi_curv(4:xpts-3)
        phi_gxx_cut = phi_gxx(4:xpts-3)
        !find average from cut arrays
        phi_curv_avg = 0.
        phi_gxx_avg = 0.
        do j=1,xpts
           phi_curv_avg = phi_curv_avg + phi_curv(j)
           phi_gxx_avg = phi_gxx_avg + phi_gxx(j)
        enddo
        phi_curv_avg = phi_curv_avg / (real(xpts))
        phi_gxx_avg = phi_gxx_avg / (real(xpts))
        
        
        ! find idex of x=0
        do j=1,xpts
           if (abs(x(j)-(0._dp))<tiny(0._dp)) then
!              print*,'x at zero=',x(j)
              zero = j
           elseif (abs(x(j)-(-120._dp))<tiny(0._dp)) then
              neg120 = j   
           elseif (abs(x(j)-(120._dp))<tiny(0._dp)) then
              pos120 = j
           elseif (abs(x(j)-(-240._dp))<tiny(0._dp)) then
              neg240 = j
           elseif (abs(x(j)-(240._dp))<tiny(0._dp)) then
              pos240 = j
           endif
!           print*,rho_x(zero)
!           delta(j) = (rho_x(j) - rho_x(zero)) / abs(rho_x(zero))
!           write(p+2,*) x(j), delta(j)
        enddo
        sum = 0.
        do j=1,xpts/2
           sum = sum + ((phi_curv(j) - phi_curv(zero)) + (phi_curv(xpts - j + 1) - phi_curv(zero)))
        enddo
        sum = sum / (xpts/2)
 
        do j=1,xpts
           delta(j) = (rho_x(j) - rho_x(zero)) / abs(rho_x(zero))
           write(p+2,*) x(j), delta(j)
        enddo
       
        phi_curv_diff = (phi_curv(neg240)-phi_curv(zero))+(phi_curv(pos240)-phi_curv(zero))
        phi_gxx_diff = (phi_gxx(neg240)-phi_gxx(zero))+(phi_gxx(pos240)-phi_gxx(zero))

        delta_phi_gxx = -(phi_gxx(pos120) - phi_gxx(zero)) / phi_gxx(zero)
        delta_phi_curv = -(phi_curv(pos120) - phi_curv(zero)) / phi_curv(zero)

        write(100,*) s, delta_phi_gxx, abs(phi_gxx(neg120)-phi_gxx(zero)), abs(phi_gxx(zero)), phi_gxx(zero), phi_curv(zero)
       
        
        write(101,*) s, phi_gxx_avg, phi_curv_avg, delta_phi_curv, delta_phi_gxx

        write(10,*) s, ham_x(zero), mom_x(zero)

        delta_rho = abs((rho_x(neg120)-rho_x(zero))) / (rho_x(zero))
!        delta_curv = abs((curv_x(neg120)-curv_x(zero))) / curv_x(zero)
        
        delta_exact = abs(-kx**2 * s**(2._dp/3._dp) / (4._dp*pi) - 2._dp * rho0)
        vel_exact = rho0 * kx * s**(-1._dp/3._dp) / (sqrt(6._dp * pi * rho0))
        curv_exact = adot * a * (1.0_dp - 2.0_dp*phi_exact) / (sqrt(1.0_dp + 2.0_dp * phi_exact))

        write(4,*) t, s, delta_rho, delta_exact, abs(vel_x(zero)), vel_exact, abs(curv_x(neg120)), curv_exact, rho_x(zero)
 
        close(p)
        close(p+1)
        close(m)
        close(m+1)
        close(m+2)
     endif
  enddo
  ! end of time-loop

  close(1)
  close(100)
  close(3)
  q=10
  do i=1,xpts
     q = q+1
     close(q)
  enddo

  deallocate(gxx_x,phi_curv,phi_gxx,x,curv_x,rho_x,delta)


     print*, 'WRITTEN: s, delta_phi_gxx, delta_phi_curv, phi_gxx(0), phi_curv(0), TO phi_growth.out'
     print*,'WRITTEN: numerical vs exact to exact_vals.out'

end program phi
