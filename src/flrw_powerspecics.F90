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
# include "fftw3.f03"

contains

    !
    ! A subroutine to generate ICs given a powerspectrum of DENSITY fluctuations; returns \phi, v^i
    !        -- for use with FLRW_Powerspectrum
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
        CCTK_REAL :: kx,ky,kz,modk2,modk_phys,modk,kspacing
        ! Real arrays for the transformed complex arrays - to be padded w/ ghosts for the output
        CCTK_REAL, dimension(:,:,:), allocatable :: delta_r,phi_r,velx_r,vely_r,velz_r

        CCTK_REAL :: Pk_interp,scalefac,Lunit
        integer :: i,j,k,n,ngh,istrt,ifin,nks,pklen
        CCTK_REAL, allocatable :: Pk(:),ks(:)
        CCTK_REAL :: C1,C3,twopi
        character(len=200) :: pkfilename

        !
        ! Define constants C1, C3 as defined (and similarly named) in eqs (23) and (25) in Macpherson et al. 2017 (in code units)
        C1 = 2.0d0 / (3.0d0 * hub**2)            ! equiv to: a_init / ( 4. * np.pi * Grhostar)
        C3 = - 2.0d0 / (3.0d0 * a0*hub)          ! equiv to: - np.sqrt( a_init / ( 6. * np.pi * Grhostar ) ) / a_init

        !   first; convert to a Fortran string check if it exists
        call CCTK_FortranString(pklen,FLRW_powerspectrum_file,pkfilename)

        ! Set up things for the power spec ICs
        call FLRW_PkICs_Setup(CCTK_ARGUMENTS,pkfilename,nx,boxlen,ks,Pk,nks,scalefac,kspacing,n,Lunit,random)

        !
        ! Loop over k-space and get k-values & interpolated power spectrum
        !    & use this to scale random field
        !
        call CCTK_INFO("Setting up Gaussian random perturbations")
        do k = 1,nx
            do j = 1,nx
                do i = 1,nx

                    ! Get our position in k-space grid
                    call FLRW_get_modk(CCTK_ARGUMENTS,i,j,k,nx,kspacing,n,kx,ky,kz,modk,modk2)
                    !
                    ! Interpolate P to this point using PHYSICAL value of modk
                    modk_phys = modk / Lunit
                    call FLRW_Interp1DLinear(modk_phys,nks,ks,Pk,Pk_interp)

                    !
                    ! Scale random field to this P(k) = delta in synchronous gauge
                    delta_sync    = random(i,j,k) * sqrt(Pk_interp/scalefac)
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


      !
      ! A subroutine to generate ICs given a spectrum of PHI perturbations; returns \phi only
      !       -- for use with FLRW_Powerspectrum_Exact
      !
      subroutine FLRW_MakePkICs_Phi(CCTK_ARGUMENTS,boxlen,nx,phi)
          implicit none
          DECLARE_CCTK_ARGUMENTS
          DECLARE_CCTK_FUNCTIONS
          DECLARE_CCTK_PARAMETERS

          CCTK_INT,  intent(in) :: nx
          CCTK_REAL, intent(in) :: boxlen
          ! Here we don't output the full array of cctk_gsh size; ghosts are added later
          CCTK_REAL, dimension(nx,nx,nx), intent(out) :: phi

          ! We make the initial data in k-space (complex arrays) and DFT back afterwards
          type(C_PTR) :: planp
          complex(C_DOUBLE_COMPLEX), dimension(nx,nx,nx) :: phi_c,random
          CCTK_REAL :: kx,ky,kz,modk2,modk_phys,modk,Lunit,kspacing
          CCTK_REAL :: Pk_interp,scalefac
          integer :: nks,i,j,k,n,pklen
          CCTK_REAL, allocatable :: Pk(:),ks(:)
          character(len=200) :: pkfilename

          !
          !   first; convert powerspectrum filename to a Fortran string to pass to setup
          call CCTK_FortranString(pklen,FLRW_phi_powerspectrum_file,pkfilename)

          !
          ! Set up things for the power spec ICs
          call FLRW_PkICs_Setup(CCTK_ARGUMENTS,pkfilename,nx,boxlen,ks,Pk,nks,scalefac,kspacing,n,Lunit,random)

          !
          ! Loop over k-space and get k-values & interpolated power spectrum
          !    & use this to scale random field
          !
          call CCTK_INFO("Setting up Gaussian random perturbations")
          do k = 1,nx
              do j = 1,nx
                  do i = 1,nx
                      !
                      ! Get k values in ~~ code units ~~
                      call FLRW_get_modk(CCTK_ARGUMENTS,i,j,k,nx,kspacing,n,kx,ky,kz,modk,modk2)

                      !
                      ! Interpolate P to this point using PHYSICAL value of modk
                      modk_phys = modk / Lunit
                      call FLRW_Interp1DLinear(modk_phys,nks,ks,Pk,Pk_interp)

                      !
                      ! Scale random field to this P(k) = phi
                      phi_c(i,j,k)  = random(i,j,k) * sqrt(Pk_interp/scalefac)
                  enddo
              enddo
          enddo
          !
          ! Transform perturbations to real space using FFTW3
          !     FFTW3 notes: https://gibbs.ccny.cuny.edu/technical/Notes/FFTW/fftw3.pdf
          !
          !  Note that FFTW3 computes the UNNORMALISED DFT
          !           --> we need to divide by nx^3 after transform
          !
          ! [ phi,phi means in-place transform (i.e. overwrites array) ]
          planp = fftw_plan_dft_3d(nx,nx,nx,phi_c,phi_c,FFTW_BACKWARD,FFTW_ESTIMATE)
          call fftw_execute_dft(planp,phi_c,phi_c)
          phi   = real(phi_c) / nx**3

          !
          ! Destroy the plan to clean up
          call fftw_destroy_plan(planp)

          deallocate(Pk,ks)

        end subroutine FLRW_MakePkICs_Phi



        !
        ! A subroutine to do the set-up to generate initial data from a given Pk
        !       -- takes the variable which stores the file name you want to read for P(k)
        !       -- reads in k,P and returns them
        !       -- sets up the 3D random field in k-space (Gaussian; unscaled)
        !       -- sets up the volume scale factor for the power spectrum
        !
        subroutine FLRW_PkICs_Setup(CCTK_ARGUMENTS,pkfilename,nx,boxlen,ks,Pk,nks,scalefac,kspacing,ni,Lunit,random)
            implicit none
            DECLARE_CCTK_ARGUMENTS
            DECLARE_CCTK_FUNCTIONS
            DECLARE_CCTK_PARAMETERS

            CCTK_INT, intent(in) :: nx         ! cubic grid spacing
            CCTK_REAL, intent(in) :: boxlen    ! cubic grid length in code units
            character(len=*), intent(in) :: pkfilename  ! string holding the path to the file name with k,P(k)

            CCTK_REAL, allocatable, intent(out) :: Pk(:),ks(:)    ! power spectrum and k-modes from pkfilename
            CCTK_REAL, intent(out) :: scalefac     ! the scale factor for P(k)
            CCTK_REAL, intent(out) :: kspacing,Lunit ! grid spacing in frequency space; length unit of box in Mpc/h
            complex(C_DOUBLE_COMPLEX), dimension(nx,nx,nx), intent(out) :: random   ! un-scaled Gaussian random field
            integer, intent(out) :: nks,ni ! number of k-values in given file; integer defining the boundary between +ve and -ve k-values

            CCTK_REAL, dimension(:,:,:), allocatable :: Rerand,Imrand
            CCTK_REAL :: dumk,dumP,spacing_phys,dV,volume
            logical :: pkexist,loop
            integer :: pkunit,ierr,nheads,i,j

            !
            ! Read the powerspectrum in from the given filename
            !
            pkexist = .False.
            inquire(file=pkfilename,exist=pkexist)
            if (pkexist) then
               ! the file exists; carry on
               open(newunit=pkunit,file=pkfilename,status='old')
            else
               ! the file does not exist; abort!!
               call CCTK_WARN(CCTK_WARN_ABORT,"Specified power spectrum file does not exist. Check your par file.")
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
            ! Allocate arrays +1 to allow for zero power at k=0 (not in the file; see doc)
            nks = nks + 1
            allocate(Pk(nks),ks(nks)); Pk = 0.d0; ks = 0.d0
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
            close(pkunit)
            !
            ! Set up scale factor to make P(k) dimensionless
            !
            volume       = FLRW_boxlength**3          ! volume of the domain in (Mpc/h)^3
            dV           = (FLRW_boxlength / nx)**3   ! volume element, i.e. dx^3, in (Mpc/h)^3
            scalefac     = dV**2/volume              ! scaling factor for Pk, in (Mpc/h)^3

            !
            ! Set up some things we need to define the k-values in the spatial loop
            kspacing   = 1.0d0 / boxlen          ! frequency spacing
            ni         = int((nx-1)/2) + 1       ! (n-1)//2 + 1: index for positive vs negative k values
            Lunit      = FLRW_boxlength / boxlen ! unit of length in Mpc/h

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
            !random = random - sum(random)/size(random)
            ! set the k=0 mode to zero consistently before scaling; replacing above
            random(1,1,1) = cmplx(0.d0,0.d0,kind(1.0d0))
            deallocate(Rerand,Imrand)

        end subroutine FLRW_PkICs_Setup



        !
        ! A subroutine to use current grid position to get the |k| for this i,j,k
        !
        subroutine FLRW_get_modk(CCTK_ARGUMENTS,i,j,k,nx,kspacing,ni,kx,ky,kz,modk,modk2)
            implicit none
            DECLARE_CCTK_ARGUMENTS
            DECLARE_CCTK_FUNCTIONS
            DECLARE_CCTK_PARAMETERS

            ! note all of this routine is in code units; scaling to Mpc/h is done outside this call

            CCTK_INT, intent(in) :: i,j,k ! current index position in the grid
            CCTK_INT, intent(in) :: nx ! the grid spacing
            CCTK_INT, intent(in) :: ni ! index separating positive vs negative k-values
            CCTK_REAL, intent(in) :: kspacing ! frequency spacing

            CCTK_REAL, intent(out) :: kx,ky,kz,modk2,modk ! k^2 and k in code units

            !
            ! Get k values in ~~ code units ~~
            if (i<=ni) then
                ! results[:N] =  arange(0, N, dtype=int)
                ! we are at the beginning of the array; as above
                kx = 2.d0 * pi * (i-1) * kspacing
            elseif (i>ni) then
                ! results[N:] = arange(-int(nx/2), 0, dtype=int)
                ! we are past the middle of the array, negative freqs
                kx = 2.d0 * pi * (- int(nx/2) + i-ni-1) * kspacing
            else
                call CCTK_WARN(CCTK_WARN_ALERT,"Invalid i value for kx")
            endif
            ! ky
            if (j<=ni) then
                ky = 2.d0 * pi * (j-1) * kspacing
            elseif (j>ni) then
                ky = 2.d0 * pi * (- int(nx/2) + j-ni-1) * kspacing
            else
                call CCTK_WARN(CCTK_WARN_ALERT,"Invalid j value for ky")
            endif
            ! kz
            if (k<=ni) then
                kz = 2.d0 * pi * (k-1) * kspacing
            elseif (k>ni) then
                kz = 2.d0 * pi * (- int(nx/2) + k-ni-1) * kspacing
            else
                call CCTK_WARN(CCTK_WARN_ALERT,"Invalid k value for kz")
            endif

            !
            ! Calculate k^2 and |k| to return
            modk2 = kx**2 + ky**2 + kz**2
            modk  = sqrt(modk2)

        end subroutine FLRW_get_modk


    end module FLRW_PowerspecICs
