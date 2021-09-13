!
! file    flrw_powerspectrum.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    A spectrum of perturbations to FLRW initial data, longitudinal gauge, zero shift
!            calls create_ics.py initial conditions generator for the given parameters
!
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine FLRW_Powerspectrum (CCTK_ARGUMENTS)
  USE init_tools
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer   :: i,j,k
  logical   :: lapse,dtlapse,shift,data,hydro
  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
  CCTK_REAL :: phi_ijk,kvalue
  !
  ! globally-size arrays (to read in initial data files)
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3))   :: phi_gs,delta_gs
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3),3) :: delta_vel_gs
  !
  ! locally-sized arrays (for this processor)
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))   :: phi, delta
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),3) :: delta_vel
  !
  CCTK_INT :: ncells(3)
  character(len=200) :: pkfilename
  character(len=200), parameter :: ierrfile = "create_ics.err" ! needs to match that set in builder.py
  integer :: dr_unit,dv_unit1,dv_unit2,dv_unit3,p_unit
  integer :: imp_ierr,pk_ierr,file_ierr,iunit,ierr
  integer :: il,jl,kl,iu,ju,ku,pkunit,pklen

  call CCTK_INFO("Initialising a POWER SPECTRUM of perturbations to an FLRW spacetime")

  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call set_logicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set parameters used in setting metric, matter parameters
  ! --> note boxlen is in code units here
  !
  call set_parameters(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)

  !
  ! Before we call the ICs generator we need to write the Pk file name to a text file
  !      --> this sucks, but it's easier than figuring out how to pass Fortran strings --> c strings --> Python strings...
  !      --> Next step: read Pk in here, pass the array into Python. Even this was a bit tricky, so I delayed it.
  !
  open(newunit=pkunit,file="pk_filename.txt",status='replace')
  ! convert CCTK string to Fortran string
  call CCTK_FortranString(pklen,FLRW_powerspectrum_file,pkfilename)
  write(pkunit,"(a)") trim(pkfilename)
  close(pkunit)

  !
  ! call initial conditions generator
  ! --> note FLRW_boxlength is in cMpc here
  !
  ! this makes the ICs and puts them into files with corresponding names as below
  !
  ! We only need to call this from one processor; save potential issues with MPI
  if (CCTK_MyProc(cctkGH)==0) then
      call CCTK_INFO("Calling create_ics for initial conditions...")
      if (ncells(1)/=ncells(2)) call CCTK_WARN(CCTK_WARN_ALERT,"non-uniform grid")

      call FLRW_SetPythonModulePath()

      call call_make_ics(a0,hub,FLRW_boxlength,boxlen(1),ncells(1),2*cctk_nghostzones(1),FLRW_random_seed)
      ! Check the error file exists, to check for errors in Python code (these don't show up nicely in CCTK output as is)
      inquire(file=ierrfile,iostat=ierr)
      if (ierr/=0) then
          call CCTK_WARN(CCTK_WARN_ALERT,"ierrorfile does not exist. Check consistency between builder.py and ierrfile. Some errors may be missing.")
      else
          ! Read in the error flags and produce relevant warnings
          !    -- order of reading here matches writing order in builder.py and create_ics.py
          open(action='read',file=ierrfile,newunit=iunit)
          read(iunit,*) imp_ierr
          ! add this warning first; because if imp_ierr/=0 then we couldn't call create_ics --> ierrfile has only one line
          if (imp_ierr/=0) call CCTK_WARN(CCTK_WARN_ABORT,"ERROR: Could not import create_ics (and/or) convert_types. ACTION: Check 'flrwsolverpath' in builder.py")
          read(iunit,*) pk_ierr
          read(iunit,*) file_ierr
          if (pk_ierr/=0)  call CCTK_WARN(CCTK_WARN_ABORT,"ERROR: Could not load powerspectrum file. ACTION: Check path in parameter file.")
          if (file_ierr/=0)  call CCTK_WARN(CCTK_WARN_ABORT,"ERROR: Could not save initial data to files.")
          close(iunit)
      endif
      call CCTK_INFO("Done making initial conditions.")
  endif
  ! All other processors wait here until we've made the ICs
  call CCTK_Barrier(ierr,cctkGH)
  
  !
  ! read in perturbations from files
  !
  open(newunit=dr_unit,file='init_delta.dat',status='old')  ! delta rho file
  open(newunit=dv_unit1,file='init_vel1.dat',status='old')  ! delta vel file [1]
  open(newunit=dv_unit2,file='init_vel2.dat',status='old')  ! delta vel file [2]
  open(newunit=dv_unit3,file='init_vel3.dat',status='old')  ! delta vel file [3]
  open(newunit=p_unit,file='init_phi.dat',status='old')     ! phi file

  !
  ! spatial loop over *global* grid size
  do k = 1, cctk_gsh(3)
     do j = 1, cctk_gsh(2)
        ! loop over ROW no. (i is COLUMN no.) (i,j,k) --> (column, row, z)

        !
        ! read perturbations from ICs files to global sized (gs) arrays
        !
        read(dr_unit,*)  delta_gs(:,j,k)
        read(dv_unit1,*) delta_vel_gs(:,j,k,1)
        read(dv_unit2,*) delta_vel_gs(:,j,k,2)
        read(dv_unit3,*) delta_vel_gs(:,j,k,3)
        read(p_unit,*)   phi_gs(:,j,k)
     enddo
  enddo
  call CCTK_INFO("Opened and read IC files")

  close(dr_unit); close(dv_unit1)
  close(dv_unit2); close(dv_unit3); close(p_unit)

  !
  ! indices for lower bound of local (processor) grid within global grid
  !
  il = cctk_lbnd(1) + 1
  jl = cctk_lbnd(2) + 1  ! indices output from cctk_lbnd start at 0, need to +1
  kl = cctk_lbnd(3) + 1
  !
  ! indices for upper bound of local grid within global grid
  !
  iu = cctk_ubnd(1) + 1
  ju = cctk_ubnd(2) + 1
  ku = cctk_ubnd(3) + 1
  !
  ! extract local part of global grid i.e. part this processor is using and store it
  !
  phi = phi_gs(il:iu, jl:ju, kl:ku)
  ! set up local grid, linear delta, vel perturbations
  delta     = delta_gs(il:iu, jl:ju, kl:ku)
  delta_vel = delta_vel_gs(il:iu, jl:ju, kl:ku, :)

  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              phi_ijk = phi(i,j,k)

              if (lapse) then
                 alp(i,j,k) = FLRW_lapse_value * sqrt(1._dp + 2._dp * phi_ijk)
              endif

              ! time deriv of lapse -- evolution of this is specified in ADMBase.
              if (dtlapse) then
                 dtalp(i,j,k) = 0._dp
              endif

              ! shift vector, always zero in this thorn
              if (shift) then
                 betax(i,j,k) = 0._dp
                 betay(i,j,k) = 0._dp
                 betaz(i,j,k) = 0._dp
              endif

              ! set perturbed metric and K_ij
              gxx(i,j,k) = asq * (1._dp - 2._dp * phi_ijk)
              gxy(i,j,k) = 0._dp
              gxz(i,j,k) = 0._dp
              gyy(i,j,k) = asq * (1._dp - 2._dp * phi_ijk)
              gyz(i,j,k) = 0._dp
              gzz(i,j,k) = asq * (1._dp - 2._dp * phi_ijk)

              kvalue     = - adot * a0 / alp(i,j,k)
              kxx(i,j,k) = kvalue * (1._dp - 2._dp * phi_ijk)
              kxy(i,j,k) = 0._dp
              kxz(i,j,k) = 0._dp
              kyy(i,j,k) = kvalue * (1._dp - 2._dp * phi_ijk)
              kyz(i,j,k) = 0._dp
              kzz(i,j,k) = kvalue * (1._dp - 2._dp * phi_ijk)

              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! perturb the matter
                 press(i,j,k) = 0._dp ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0._dp
                 ! set up linearly-perturbed matter data
                 rho(i,j,k)   = rho0 * (1._dp + delta(i,j,k))
                 vel(i,j,k,:) = delta_vel(i,j,k,:)

              endif

           endif

        enddo
     enddo
  enddo

  !
  ! make sure the metric type is physical (not conformal)
  !
  call check_metric()


end subroutine FLRW_Powerspectrum
