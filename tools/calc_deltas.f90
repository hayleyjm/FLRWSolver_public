!program to calculate ONLY deltaRHO, using file-splitter output files
program delta_rho
  implicit none

  integer,parameter:: n=int(1e5),dp=8
  real(dp),allocatable,dimension(:):: rho, x, vel, phi, gxx, curv
  real(dp):: deltarho, t, dum, maxrho, minrho, avgrho, deltav, maxv, minv, deltaphi, maxphi, minphi
  real(dp):: deltacurv, maxcurv, mincurv, avgcurv, avgphi
  real(dp) :: tf, tl(n), alp(n), sl, sf, dtl, tl_temp, alp_temp, a, asq_flrw, asq_temp
  real(dp), parameter:: pi=3.141592653589793238462643383279502884197, rho0=1.e-8, a0=1.
  real(dp) :: gxx_max, gxx_min, gxxmax_temp, gxxmin_temp
  integer:: i, number, res, npts, m, it, dit, itmin, header, j, l, count, itcac, lapse, p
  character(len=40):: file_rho, file_v, file_gxx, file_curv
!  character(len=100):: lapse
  logical:: reading, nopc
!  dit = 512
  itmin = 0
  

  print*,'Look at every nth iteration, enter n:'
  read(*,*) number

  print*,'Enter resolution of run (dn):'
  read(*,*) res
  npts = 480 / res ! number of grid points based on resolution

  allocate(rho(npts),x(npts),vel(npts),phi(npts),gxx(npts),curv(npts))

  print*,'Are there headers in your files? 1=yes, 2=no'
  read(*,*) header

  m = 20

! check whether files are in PostCactus format or not
  inquire(file='admbase::lapse.maximum.asc',exist=nopc) ! if this exists, then nopc=.True.
  if (nopc) then
     ! open old format lapse file for converting tl --> tf
     open(unit=12,file='admbase::lapse.maximum.asc',status='old')
     open(unit=13,file='admbase::metric.maximum.asc',status='old')
     open(unit=14,file='admbase::metric.minimum.asc',status='old')
     ! use delta calculation the old way
     dit = 512 ! old dit
  else 
     ! open Postcactus format file
     open(unit=12,file='alp.maximum.asc',status='old')
     open(unit=13,file='gxx.maximum.asc',status='old')! for reconstructing phimax, phimin
     open(unit=14,file='gxx.minimum.asc',status='old')
     open(unit=15,file='dens.maximum.asc',status='old') ! for calculating delta rho
     open(unit=16,file='dens.minimum.asc',status='old')
     dit = 1 ! new dit
  endif


  ! open FRIEDMANN gxx file for reconstructing phi, not considering finite diff errors
  ! open(unit=13,file='../FRIEDMANN_GXX_LAPSE.asc',status='old')
!  print*, 'Enter gauge used: harmonic(1), 1+log(2) or 1+log-old(3)'
!  read(*,*) lapse

!  if (lapse==1) then !harmonic gauge
!     open(unit=13,file='../FRIEDMANN_GXX_HARMONIC.asc',status='old')
!  elseif (lapse==2) then ! 1+log with f=2/alp
!     open(unit=13,file='../FRIEDMANN_GXX_1PLOG.asc',status='old')
!  elseif (lapse==3) then ! 1+log with f=1/alp
!     open(unit=13,file='../FRIEDMANN_GXX_1PLOGOLD.asc',status='old')
!  endif

! ^^ THIS CAN'T WORK BECAUSE FLRW FILES HAVE DTFAC=10... EXACT A USED INSTEAD


! read out headers of .max files
  if (nopc) then
     p = 10
  else
     p = 9
  endif
  do i=1,p
     read(12,*) ! 10 header lines in old format files, 9 in new format
     read(13,*)
     read(14,*)
  enddo

! open file for writing delta rho
  open(unit=10,file='deltas.out',status='replace')
  write(10,*) '#tf  sf   deltarho   deltav  deltaphi   avgphi'

  alp(1) = 1.0 ! initial lapse

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
        read(13,*) dum, dum, gxxmax_temp
        read(14,*) dum, dum, gxxmin_temp
        if (itcac==it) then
           tl(i) = tl_temp
           alp(i) = alp_temp
           gxx_max = gxxmax_temp
           gxx_min = gxxmin_temp
           reading = .False.
        endif
     enddo

! calculate friedmann time from lapse time
     if (i/=1) then
        dtl = tl(i) - tl(i-1)
        tf = tf + 0.5_dp * dtl * (alp(i-1) + alp(i))
     endif
! calculate 's' and FLRW scale factor a
     sf = 1._dp + sqrt(6._dp * pi * rho0) * tf
     a = a0 * sf**(2._dp/3._dp) ! exact a, if commented out: using FLRW a as above
     asq_flrw = a**2

! open existing density(x) files to calculate deltarho    
     write(file_rho,'(a,i9.9,a)')'density_',it,'.dat'
     open(unit=m,file=file_rho,status='old')
! open existing vel(x) files to calculate deltav
     write(file_v,'(a,i9.9,a)')'velocity_',it,'.dat'
     open(unit=m+1,file=file_v,status='old')
! open existing gxx(x) files to reconstruct phi
     write(file_gxx,'(a,i9.9,a)')'metric_',it,'.dat'
     open(unit=m+2,file=file_gxx,status='old')
! open existing curv(x) files 
!     write(file_curv,'(a,i9.9,a)')'curv_',it,'.dat'
!     open(unit=m+3,file=file_curv,status='old')

     if (header==1) then
        do j=1,6 ! read out headers from regular file splitter
           read(m,*)
           read(m+1,*)
           read(m+2,*)
 !          read(m+3,*)
        enddo
     else ! read out timestamp line on files from mpi file splitter
        read(m,*)
        read(m+1,*)
        read(m+2,*)
  !      read(m+3,*)
     endif

! loop over x values and read x, rho, v, gxx into arrays & calculate phi
     do j=1,npts
        read(m,*) (dum,l=1,8), t, x(j), dum, dum, rho(j)
        read(m+1,*) (dum,l=1,12), vel(j)
        read(m+2,*) (dum,l=1,12), gxx(j)
   !     read(m+3,*) (dum,l=1,12), curv(j)
        phi(j) = 0.5_dp * (1._dp - gxx(j) / asq_flrw) ! calculate phi(x) from gxx.x slice along y=z=0 (this ISN'T maximum of phi - this is in centre of domain with 3d perturbation!!)
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

! calculate deltav
     deltav = 0.5_dp * (maxv - minv)

! calculate phi_max and min
     maxphi = maxval(phi) ! this max, mins & avg for use with phi(x) from 1D gxx slice
     minphi = minval(phi)
     avgphi = sum(phi) / size(phi)

!     maxphi = 0.5_dp * (1._dp - gxx_max / asq_flrw) ! these max, mins & avg for use from gxx_max and min
!     minphi = 0.5_dp * (1._dp - gxx_min / asq_flrw)
!     avgphi = maxphi - deltaphi

     deltaphi = 0.5_dp * (maxphi - minphi)

! calculate curvature
 !    maxcurv = maxval(curv)
 !    mincurv = minval(curv)
 !    avgcurv = sum(curv) / size(curv)
 !    deltacurv = 0.5_dp * (maxcurv - mincurv) / avgcurv

     print*,'WRITING tl=',tl(i),'tf=',tf,' deltarho=',deltarho, 'deltaphi=',abs(deltaphi)
     write(10,*) tf, sf, deltarho, deltav,  abs(deltaphi), abs(avgphi)

     close(m)
     close(m+1)
     close(m+2)
     close(m+3)
  enddo

close(12)
close(13)
close(14)
end program delta_rho
