!program to calculate ONLY delta_RHO and delta_VEL, and convert time from lapse time for exact solution comparison. writes both times and corresponding 's' values, deltarho and deltav to output file
!! REQUIRES: --> lapse.max, density_*, velocity_*

program delta_rho
  implicit none

  integer,parameter:: n=int(1e5),dp=8
  real(dp),allocatable,dimension(:):: rho, x, vel
  real(dp):: deltarho, t, dum, maxrho, minrho, avgrho, deltav, maxv, minv
  real(dp) :: tf, tl(n), alp(n), sl, sf, dtl, tl_temp, alp_temp
  real(dp), parameter:: pi=3.141592653589793238462643383279502884197, rho0=1.e-8
  integer:: i, number, res, npts, m, it, dit, itmin, header, j, l, count, itcac
  character(len=40):: file_rho, file_v
  logical:: reading
  dit = 512
  itmin = 0

  print*,'Look at every nth iteration, enter n:'
  read(*,*) number

  print*,'Enter resolution of run (dn):'
  read(*,*) res
  npts = 480 / res ! number of grid points based on resolution

  allocate(rho(npts),x(npts),vel(npts))

  print*,'Are there headers in your files? 1=yes, 2=no'
  read(*,*) header

  m = 20

! open lapse file for converting tl --> tf
  open(unit=12,file='admbase::lapse.maximum.asc',status='old')
  do i=1,10
     read(12,*)
  enddo

! open file for writing delta rho
  open(unit=10,file='delta_rho.out',status='replace')
  write(10,*) '# tl              tf            sl             sf             deltarho           deltav'

  tf = 0.
  count = 0
  do i=1,n
     m = m+1
     count = count + 1

     ! calculate iteration
     if (i==1) then
        it = itmin
     else
        it = itmin + ((i-1)*dit*number - dit)
     endif
 
! read lapse time, alpha every number^th iteration
     reading = .True.
     do while (reading)
        read(12,*) itcac, tl_temp, alp_temp
        if (itcac==it) then
           tl(i) = tl_temp
           alp(i) = alp_temp
           reading = .False.
        endif
     enddo

! calculate friedmann time from lapse time
     if (i/=1) then
        dtl = tl(i) - tl(i-1)
        tf = tf + 0.5_dp * dtl * (alp(i-1) + alp(i))
     endif

! open existing density(x) files to calculate deltarho    
     write(file_rho,'(a,i9.9,a)')'density_',it,'.dat'
     open(unit=m,file=file_rho,status='old')
! open existing vel(x) files to calculate deltav
     write(file_v,'(a,i9.9,a)')'velocity_',it,'.dat'
     open(unit=m+1,file=file_v,status='old')

     if (header==1) then
        do j=1,6 ! read out headers from regular file splitter
           read(m,*)
           read(m+1,*)
        enddo
     else ! read out timestamp line on files from mpi file splitter
        read(m,*)
        read(m+1,*)
     endif

! loop over x values and read x, rho into arrays
     do j=1,npts
        read(m,*) (dum,l=1,8), t, x(j), dum, dum, rho(j)
        read(m+1,*) (dum,l=1,12), vel(j)
        if (x(j)==0) then
           avgrho = rho(j)
        elseif (x(j)==-120) then
           maxrho = rho(j)
        endif
     enddo
     if (t/=tl(i)) then
        print*,'ERROR: times from density files dont match time from alp file !!!!!'
        print*,t, tl(i)
     endif

! calculate max, min, avg of density array
     maxrho = maxval(rho)
     minrho = minval(rho)
     avgrho = sum(rho) / size(rho)

! calculate max, min of vel
     maxv = maxval(vel)
     minv = minval(vel)

! calculate deltarho
     deltarho = 0.5_dp * (maxrho - minrho) / avgrho
!     deltarho = (maxrho - avgrho) / avgrho

! calculate deltav
     deltav = 0.5_dp * (maxv - minv)

     sl = 1._dp + sqrt(6._dp * rho0 * pi) * t ! "lapse s"
     sf = 1._dp + sqrt(6._dp * rho0 * pi) * tf ! "friedmann s"

     write(10,*) t, tf, sl, sf, deltarho, deltav

     close(m)
     close(m+1)
  enddo

end program delta_rho
