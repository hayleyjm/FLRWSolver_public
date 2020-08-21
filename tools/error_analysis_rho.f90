!! program to calculate the error in delta_rho at each timestep, using exact solution. also calculates error in rho(x) at each timestep and writes to individual error_density_* files.
!! REQUIRES: --> delta_vals (as calculated from calc_deltas.f90), density_*
program error_analysis_rho
  implicit none

  integer, parameter:: n=int(1e5),dp=8
  real(dp), parameter:: pi=3.141592653589793238462643383279502884197, rho0=1.e-8, a0=1., kx=2.*pi/480.
  real(dp):: delta, s, delta_exact, dum, delta_error, s2, t, delta_exact2, a, rhobar, maxerror
  real(dp), allocatable, dimension(:):: rho_x, error_x, rho_exact_x, x
  integer:: i, it, itmin, npts, number, res, m, j, l, dit, header
  character(len=40):: file_rho, file_error, file_rho_exact

  itmin = 0
  dit = 512
  m = 50

  print*,'Look at every nth iteration, enter n:'
  read(*,*) number

  print*,'Enter resolution of run:'
  read(*,*) res
  npts = 480 / res ! number of grid points based on resolution

  allocate(error_x(npts),rho_x(npts),rho_exact_x(npts),x(npts))

  print*,'Are there headers in your files? 1=yes, 2=no'
  read(*,*) header

!open file to read delta rho
  open(unit=10,file='delta_vals.out',status='old')
  read(10,*) ! read out header line 

!open file to write delta(t) error to
  open(unit=12,file='error_delta.dat',status='replace')
  write(12,*) '# s      delta       delta_exact           error'

!open file to write max error of rho(x) over time
  open(unit=13,file='error_max_rhox.dat',status='replace')
  write(13,*) '# s         max_error of rho_x'


  do i=1,n+1

     m = m+1

     ! calculate iteration
     if (i==1) then
        it = itmin
     else
        it = itmin + ((i-1)*dit*number - dit)
     endif
     print*,'ITERATION is', it
     
! read cactus delta and calculate exact value, find error and write to file
     read(10,*) s, delta
     delta_exact = - kx**2 / (4._dp * pi) * s**(2._dp/3._dp) - 2._dp * rho0     
     delta_error = (abs(delta) - abs(delta_exact)) / abs(delta_exact)
     write(12,*) s, delta, abs(delta_exact), abs(delta_error)


!open existing rho(x) files
     write(file_rho,'(a,i9.9,a)')'density_',it,'.dat'
     open(unit=m,file=file_rho,status='old')
     if (header==1) then
        do j=1,6 ! headers from original (non-mpi) file splitter
           read(m,*)
        enddo
     else ! headers from mpi file splitter, timestamp row only
        read(m,*)
     endif

!write files for exact density
!     write(file_rho_exact,'(a,i9.9,a)')'exact_density_',it,'.dat'
!     open(unit=m+5,file=file_rho_exact,status='replace')

!create files to write error(x) 
     write(file_error,'(a,i9.9,a)')'error_density_',it,'.dat'
     open(unit=m+2,file=file_error,status='replace')

! calculate error at each point in x for rho
     do j=1,npts
        read(m,*) (dum,l=1,8), t, x(j), dum, dum, rho_x(j)
        s2 = 1._dp + sqrt(6._dp * pi * rho0) * t
        delta_exact2 = - kx**2 / (4._dp * pi) * s2**(2._dp/3._dp) - 2._dp * rho0
        
        a = a0 * s2**(2._dp/3._dp)
        rhobar = rho0 * a0**3 / (a**3)
        rho_exact_x(j) = rhobar * (1._dp + delta_exact * sin(kx *x(j)))

        error_x(j) = (rho_x(j) - rho_exact_x(j)) / rho_exact_x(j)

        write(m+2,*) x(j), error_x(j)
 !       write(m+5,*) x(j), rho_exact_x(j)
     enddo

     maxerror = maxval(error_x)
     write(13,*) s2, maxerror

     close(m)
     close(m+2)
  enddo

close(12)
close(13)
end program error_analysis_rho
