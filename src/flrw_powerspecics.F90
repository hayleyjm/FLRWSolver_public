! file    powerspec_ics.F90
! author  Hayley Macpherson
! date    05.08.2022
! desc    A module to generate Gaussian random initial data from a given powerspectrum
!
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module FLRW_PowerspecICs
  use, intrinsic :: iso_c_binding ! we need this for the fftw calls
  use FLRW_InitTools, only: pi,FLRW_Interp1DLinear,FLRW_GetRandomNormal3D
  implicit none
  include 'fftw3.f03'

contains

    !
    ! A subroutine to generate the ICs
    !
    subroutine FLRW_MakePkICs(CCTK_ARGUMENTS,a0,hub,boxlen,nx,delta_gs,phi_gs,velx_gs,vely_gs,velz_gs)
        implicit none
        DECLARE_CCTK_ARGUMENTS
        DECLARE_CCTK_FUNCTIONS
        DECLARE_CCTK_PARAMETERS

        CCTK_INT,  intent(in) :: nx
        CCTK_REAL, intent(in) :: a0,hub,boxlen
        ! Global size grid arrays to output (including boundaries)
        CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)), intent(out) :: delta_gs,phi_gs,velx_gs,vely_gs,velz_gs

        ! We make the initial data in k-space (complex arrays) and DFT back afterwards
        type(C_PTR) :: pland,planp,planvx,planvy,planvz
        complex(C_DOUBLE_COMPLEX), dimension(nx,nx,nx) :: delta,phi,velx,vely,velz,random
        complex(C_DOUBLE_COMPLEX) :: delta_sync
        CCTK_REAL :: kx,ky,kz,modk2,modk_phys,modk,Lunit,val
        ! Real arrays for the transformed complex arrays - to be padded w/ ghosts for the output
        CCTK_REAL, dimension(:,:,:), allocatable :: delta_r,phi_r,velx_r,vely_r,velz_r
        CCTK_REAL, dimension(:,:,:), allocatable :: Rerand,Imrand

        CCTK_REAL :: Pk_interp,Pkscale,volume,dV,scale
        integer :: pklen,nks,nheads,ierr,i,j,k,n,pkunit,ngh,istrt,ifin
        CCTK_REAL, allocatable :: Pk(:),ks(:)
        CCTK_REAL :: dumk,dumP,spacing,spacing_phys,C1,C3,twopi
        logical :: loop,pkexist
        character(len=200) :: pkfilename        

        !
        ! Define constants C1, C3 from Macpherson et al. 2017 (in code units)
        C1 = 2.0d0 / (3.0d0 * hub**2)            ! equiv to: a_init / ( 4. * np.pi * Grhostar)
        C3 = - 2.0d0 / (3.0d0 * a0*hub)          ! equiv to: - np.sqrt( a_init / ( 6. * np.pi * Grhostar ) ) / a_init

        !
        ! Read the powerspectrum in from the given filename
        !
        !
        !   first; convert to a Fortran string check if it exists
        call CCTK_FortranString(pklen,FLRW_powerspectrum_file,pkfilename)
        pkexist = .False.
        inquire(file=pkfilename,exist=pkexist)
        if (pkexist) then
           ! the file exists; carry on
           open(newunit=pkunit,file=pkfilename,status='old')
        else
           ! the file does not exist; abort!!
           call CCTK_WARN(CCTK_WARN_ABORT,"The power spectrum file you specified does not exist. Check FLRWSolver::FLRW_powerspectrum_file and use the example file in FLRWSolver/powerspectrum/ or generate your own using CLASS.")
        endif
        !
        ! We need to count the number of lines in the file
        !     (typically ~ 100 lines, this will be quick)
        ierr = 0; nks = 0; nheads = 0
        loop = .True.
        do while(loop) ! skip the first instances of ierr/=0 due to header
            read(pkunit,*,iostat=ierr) dumk, dumP
            if (ierr==0) then
                nks = nks + 1 ! if real number ; add to count, else; header/end of file
            else
                if (nks>0) then
                    ! ierr /=0 and we've counted k's --> end of file
                    loop = .False.
                else
                    ! ierr /=0 and nks==0 --> still at header, count them
                    nheads = nheads + 1
                endif
            endif
        enddo
        ! Allocate arrays +1 to allow for zero power at k=0
        allocate(Pk(nks+1),ks(nks+1)); Pk = 0.d0; ks = 0.d0
        !
        ! Now read in the data to these arrays
        rewind(pkunit) ! rewind to the beginning of the file
        call CCTK_INFO("Reading power spectrum")
        j = 2 ! start at 2 to keep Pk(1) = ks(1) = 0
        do i=1,nks+nheads
            read(pkunit,*,iostat=ierr) dumk, dumP
            if (ierr==0) then
                ! read was successful; store
                ks(j) = dumk
                Pk(j) = dumP
                j = j + 1
            !else; read was not successful; must be header & do nothing
            endif
        enddo
        !
        ! Set up some things about the domain we need for the k-arrays, Pk scaling, etc
        Lunit   = FLRW_boxlength / boxlen  ! unit of length
        spacing = boxlen / nx              ! dx in code units
        spacing_phys = spacing * Lunit     ! dx in physical units (Mpc/h)

        !
        ! Set up complex Gaussian random array
        !      and scale to be centered around zero (this will be delta)
        !
        allocate(Rerand(nx,nx,nx),Imrand(nx,nx,nx))
        ! Note that two calls to this generator with same seed will generate the same numbers
        !       so we change the seed slightly, which will still be consistent between runs
        call FLRW_GetRandomNormal3D(nx,nx,nx,Rerand,FLRW_random_seed)
        call FLRW_GetRandomNormal3D(nx,nx,nx,Imrand,FLRW_random_seed*42+1068)
        random = cmplx(Rerand,Imrand,kind(1.0d0))
        random = random - sum(random)/size(random)
        deallocate(Rerand,Imrand)
        !
        ! Set up scale factor to make P(k) dimensionless
        !
        volume = FLRW_boxlength**3   ! volume of the domain in (Mpc/h)^3
        dV     = spacing_phys**3     ! volume element, i.e. dx^3, in (Mpc/h)^3
        scale  = dV**2/volume        ! scaling factor for Pk, in (Mpc/h)^3

        !
        ! Loop over k-space and get k-values & interpolated power spectrum
        !    & use this to scale random field
        !

        ! Set up some things we need to define the k-values
        val   = 1.0d0 / (nx * spacing)  ! frequency spacing
        n     = int((nx-1)/2) + 1       ! (n-1)//2 + 1: index for positive vs negative k values
        twopi = 2.0d0 * pi
        call CCTK_INFO("Setting up Gaussian random perturbations")
        do k = 1,nx
            do j = 1,nx
                do i = 1,nx
                    !
                    ! Get k values in ~~ code units ~~
                    if (i<=n) then
                        ! results[:N] =  arange(0, N, dtype=int)
                        ! we are at the beginning of the array; as above
                        kx = twopi * (i-1) * val
                    elseif (i>n) then
                        ! results[N:] = arange(-int(nx/2), 0, dtype=int)
                        ! we are past the middle of the array, negative freqs
                        kx = twopi * (- int(nx/2) + i-n-1) * val
                    else
                        call CCTK_WARN(CCTK_WARN_ALERT,"Invalid i value for kx")
                    endif
                    ! ky
                    if (j<=n) then
                        ky = twopi * (j-1) * val
                    elseif (j>n) then
                        ky = twopi * (- int(nx/2) + j-n-1) * val
                    else
                        call CCTK_WARN(CCTK_WARN_ALERT,"Invalid j value for ky")
                    endif
                    ! kz
                    if (k<=n) then
                        kz = twopi * (k-1) * val
                    elseif (k>n) then
                        kz = twopi * (- int(nx/2) + k-n-1) * val
                    else
                        call CCTK_WARN(CCTK_WARN_ALERT,"Invalid k value for kz")
                    endif
                    ! Calculate k^2 and |k|
                    modk2 = kx**2 + ky**2 + kz**2
                    modk  = sqrt(modk2)
                    !
                    ! Interpolate P to this point using PHYSICAL value of modk
                    modk_phys = modk / Lunit
                    call FLRW_Interp1DLinear(modk_phys,nks,ks,Pk,Pk_interp)

                    !
                    ! Scale random field to this P(k) = delta in synchronous gauge
                    Pkscale       = sqrt(Pk_interp/scale)
                    delta_sync    = random(i,j,k) * Pkscale
                    !
                    ! Calculate phi from delta in synchronous
                    if (modk2==0.0d0) then
                        ! k=0 point; set to zero to avoid NaN/infs
                        phi(i,j,k) = cmplx(0.0d0,0.0d0,kind(1.0d0))
                    else
                        ! k/=0 points
                        phi(i,j,k) = - delta_sync / (C1 * modk2)
                    endif
                    !
                    ! Calculate delta in longitudinal gauge & delta_vel from phi
                    !    (all in code units)
                    delta(i,j,k)  = - phi(i,j,k) * (C1 * modk2 + 2.0d0)
                    velx(i,j,k)   = cmplx(0.0d0,C3*kx,kind(1.0d0)) * phi(i,j,k) ! C3 * 1j * kx * phi
                    vely(i,j,k)   = cmplx(0.0d0,C3*ky,kind(1.0d0)) * phi(i,j,k) ! C3 * 1j * ky * phi
                    velz(i,j,k)   = cmplx(0.0d0,C3*kz,kind(1.0d0)) * phi(i,j,k) ! C3 * 1j * kz * phi
                enddo
            enddo
        enddo
        !
        ! Transform perturbations to real space using FFTW3
        !     FFTW3 notes: https://gibbs.ccny.cuny.edu/technical/Notes/FFTW/fftw3.pdf
        !
        ! Allocate real-space arrays
        allocate(delta_r(nx,nx,nx),phi_r(nx,nx,nx),velx_r(nx,nx,nx),&
            & vely_r(nx,nx,nx),velz_r(nx,nx,nx))
        !
        !  Note that FFTW3 computes the UNNORMALISED DFT
        !           --> we need to divide by nx^3 after transform
        !
        ! [ phi,phi means in-place transform (i.e. overwrites array) ]
        planp = fftw_plan_dft_3d(nx,nx,nx,phi,phi,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(planp,phi,phi)
        phi_r = real(phi) / nx**3

        pland = fftw_plan_dft_3d(nx,nx,nx,delta,delta,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pland,delta,delta)
        delta_r = real(delta) / nx**3

        planvx = fftw_plan_dft_3d(nx,nx,nx,velx,velx,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(planvx,velx,velx)
        velx_r = real(velx) / nx**3

        planvy = fftw_plan_dft_3d(nx,nx,nx,vely,vely,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(planvy,vely,vely)
        vely_r = real(vely) / nx**3

        planvz = fftw_plan_dft_3d(nx,nx,nx,velz,velz,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(planvz,velz,velz)
        velz_r = real(velz) / nx**3
        !
        ! Destroy all the plans to clean up
        call fftw_destroy_plan(pland)
        call fftw_destroy_plan(planp)
        call fftw_destroy_plan(planvx)
        call fftw_destroy_plan(planvy)
        call fftw_destroy_plan(planvz)

        !
        ! Put data into global sized grid arrays (including boundary zones) to be passed out
        !
        ! Initialise to zero such that boundary zones remain zero
        delta_gs = 0.0d0; phi_gs = 0.0d0; velx_gs = 0.0d0; vely_gs = 0.0d0; velz_gs = 0.0d0
        ! Fill interior values
        !   Define indexing to extract interior zone
        ngh   = cctk_nghostzones(1) ! uniform grid only
        istrt = ngh+1; ifin = nx+ngh
        !  Set interior zone
        delta_gs(istrt:ifin,istrt:ifin,istrt:ifin) = delta_r
        phi_gs(istrt:ifin,istrt:ifin,istrt:ifin)   = phi_r
        velx_gs(istrt:ifin,istrt:ifin,istrt:ifin)  = velx_r
        vely_gs(istrt:ifin,istrt:ifin,istrt:ifin)  = vely_r
        velz_gs(istrt:ifin,istrt:ifin,istrt:ifin)  = velz_r

        deallocate(delta_r,phi_r,velx_r,vely_r,velz_r,Pk,ks)

      end subroutine FLRW_MakePkICs


    end module FLRW_PowerspecICs
