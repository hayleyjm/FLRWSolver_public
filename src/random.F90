module random
  !
  ! A module to do some things related to generating random fields
  !
  use iso_fortran_env, only: int64
  use init_tools, only:dp,pi,c_double
  implicit none

contains

    !
    ! A subroutine to draw samples from a normal distribution
    !     using the Box-Muller transformation from a distribution in [0,1]
    !
    subroutine get_random_normal_3d(nx,ny,nz,randnums,rseed)
        integer, intent(in) :: nx,ny,nz
        real(c_double), intent(out) :: randnums(nx,ny,nz)
        integer, intent(in), optional :: rseed
        real(c_double) :: rand1(nx,ny,nz),rand2(nx,ny,nz)

        ! 1. Generate the random numbers in [0,1]
        if (present(rseed)) then
            call get_randnums_3d(nx,ny,nz,rand1,rseed)
            call get_randnums_3d(nx,ny,nz,rand2,rseed)
        else
            call get_randnums_3d(nx,ny,nz,rand1)
            call get_randnums_3d(nx,ny,nz,rand2)
        endif
        ! 2. Now get the normally distributed numbers from these (see Wikipedia page)
        randnums = sqrt(-2._dp*log(rand1))*cos(2._dp*pi*rand2)

    end subroutine get_random_normal_3d


  !
  ! A subroutine to get randomly-drawn 3D array of numbers
  !
  subroutine get_randnums_3d(nx,ny,nz,randnums,rseed)
      integer, intent(in) :: nx,ny,nz
      real(c_double), intent(out) :: randnums(nx,ny,nz)
      integer, intent(in), optional :: rseed

      if (present(rseed)) then
          call init_random_seed(rseed)
      else
          call init_random_seed()
      endif
      call RANDOM_NUMBER(randnums)

  end subroutine get_randnums_3d



  !
  ! A subroutine to initialise the random_number generator
  !    picked from here: https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran/RANDOM_005fSEED.html
  !
  subroutine init_random_seed(rseed)
    integer, intent(in), optional :: rseed ! chosen random seed for repeatability
    integer, allocatable :: seed(:)
    integer :: i,n,un,istat,dt(8),pid
    integer(int64) :: t,rs

    call random_seed(size=n)
    allocate(seed(n))
    if (present(rseed)) then
      !
      ! We've got an rseed, use that
      !
      !seed = rseed
      rs = rseed
      do i=1,n
         seed(i) = PRNG_simp(rs)
      enddo
    else
      ! we haven't been given a seed; find one
      !
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=istat)
      if (istat==0) then
         read(un) seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         call system_clock(t)
         if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                 + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                 + dt(3) * 24_int64 * 60 * 60 * 1000 &
                 + dt(5) * 60 * 60 * 1000 &
                 + dt(6) * 60 * 1000 + dt(7) * 1000 &
                 + dt(8)
         endif
         pid = getpid()
         t   = ieor(t,int(pid,kind(t)))
         do i=1,n
            seed(i) = PRNG_simp(t)
         enddo
      endif
    endif
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine init_random_seed



  !
  ! A simple PRNG found on stackexchange: https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran/RANDOM_005fSEED.html
  !  --> gives a PRN (pseudo-random number) given a single integer...
  !  --> SOME repeated numbers, e.g. calling with seed for the 10th time
  !       will give the same number as calling with seed*10 for the 1st time
  !
  integer function PRNG_simp(s)
    use iso_fortran_env, only: int64
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    PRNG_simp = int(mod(s, int(huge(0), int64)), kind(0))
  end function PRNG_simp



end module random
