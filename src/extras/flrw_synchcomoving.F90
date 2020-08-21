! file    flrw_synchcomoving.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    A single mode perturbation (for now) in the synchronous (alp=1), comoving gauge. See Macpherson & Bruni overleaf draft
!               note: these ICs are solved from the constraints directly (not linear perturb. theory) so they are generated externally, and read in from files.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine FLRW_SynchComoving (CCTK_ARGUMENTS)
  USE init_tools
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer   :: i,j,k
  logical   :: lapse,dtlapse,shift,data,hydro
  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen(3)
  CCTK_REAL :: kx,ky,kz,modk
  !
  ! globally-size arrays (to read in initial data files)
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: delta_gs,rc_gs!,chi_gs
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)) :: dxdxchi_gs,dxdychi_gs,dxdzchi_gs,dydychi_gs,dydzchi_gs,dzdzchi_gs
  !
  ! locally-sized arrays (for this processor)
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: delta,rc!,chi
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: dxdxchi,dxdychi,dxdzchi,dydychi,dydzchi,dzdzchi
  CCTK_INT  :: ncells(3)
  !
  character(len=100) :: deltafile,rcfile,ics_dir!,chifile
  character(len=100) :: dxdxchifile,dxdychifile,dxdzchifile,dydychifile,dydzchifile,dzdzchifile
  integer :: dr_unit,rc_unit,xx_unit,xy_unit,xz_unit,yy_unit,yz_unit,zz_unit!,chi_unit
  integer :: ics_dir_len,il,jl,kl,iu,ju,ku

  call CCTK_INFO("Initialising a perturbation to an FLRW spacetime in the SYNCHRONOUS, COMOVING gauge")
  
  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call set_logicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set parameters used in setting metric, matter parameters
  !
  call set_parameters(CCTK_ARGUMENTS,a0,rho0,asq,rhostar,hub,adot,hubdot,boxlen,ncells)

  ! wavenumbers for each direction (all the same since we assume a regular grid)
  kx = 2._dp * pi / boxlen(1)
  ky = 2._dp * pi / boxlen(2)
  kz = 2._dp * pi / boxlen(3)
  modk = sqrt(kx**2 + ky**2 + kz**2)
  
  !
  ! open files to read delta, chi, rc perturbs from files
  !

  ! convert CCTK_STRING "FLRW_ICs_dir" to Fortran string
  call CCTK_FortranString(ics_dir_len,FLRW_ICs_dir,ics_dir)

  write(deltafile,'(a,a)')trim(ics_dir),'/delta.dat'
  !write(chifile,'(a,a)')trim(ics_dir),'/chi.dat'
  write(rcfile,'(a,a)')trim(ics_dir),'/rc.dat'
  write(dxdxchifile,'(a,a)')trim(ics_dir),'/dxdxchi.dat'
  write(dxdychifile,'(a,a)')trim(ics_dir),'/dxdychi.dat'
  write(dxdzchifile,'(a,a)')trim(ics_dir),'/dxdzchi.dat'
  write(dydychifile,'(a,a)')trim(ics_dir),'/dydychi.dat'
  write(dydzchifile,'(a,a)')trim(ics_dir),'/dydzchi.dat'
  write(dzdzchifile,'(a,a)')trim(ics_dir),'/dzdzchi.dat'

  open(newunit=dr_unit,file=deltafile,status='old')   ! delta rho file
  !open(newunit=chi_unit,file=chifile,status='old')    ! chi file file [1]
  open(newunit=rc_unit,file=rcfile,status='old')      ! R_c file [2]
  open(newunit=xx_unit,file=dxdxchifile,status='old') ! d2dx2(chi) file
  open(newunit=xy_unit,file=dxdychifile,status='old') ! d2dxdy(chi) file
  open(newunit=xz_unit,file=dxdzchifile,status='old') ! d2dxdx(chi) file
  open(newunit=yy_unit,file=dydychifile,status='old') ! d2dy2(chi) file
  open(newunit=yz_unit,file=dydzchifile,status='old') ! d2dydz(chi) file
  open(newunit=zz_unit,file=dzdzchifile,status='old') ! d2dz2(chi) file

  do k = 1, cctk_gsh(3)
     do j = 1, cctk_gsh(2) 
        ! loop over ROW no. (i is COLUMN no.) (i,j,k) --> (column, row, z)

        !
        ! read perturbations from 3D files to global sized (gs) arrays
        !
        read(dr_unit,*)  delta_gs(:,j,k)
        !read(chi_unit,*) chi_gs(:,j,k)
        read(rc_unit,*)  rc_gs(:,j,k)
        read(xx_unit,*)  dxdxchi_gs(:,j,k)
        read(xy_unit,*)  dxdychi_gs(:,j,k)
        read(xz_unit,*)  dxdzchi_gs(:,j,k)
        read(yy_unit,*)  dydychi_gs(:,j,k)
        read(yz_unit,*)  dydzchi_gs(:,j,k)
        read(zz_unit,*)  dzdzchi_gs(:,j,k)
     enddo
  enddo
  call CCTK_INFO("Opened and read IC files")

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
  delta   = delta_gs(il:iu, jl:ju, kl:ku)
  !chi     = chi_gs(il:iu, jl:ju, kl:ku)
  rc      = rc_gs(il:iu, jl:ju, kl:ku)
  dxdxchi = dxdxchi_gs(il:iu, jl:ju, kl:ku)
  dxdychi = dxdychi_gs(il:iu, jl:ju, kl:ku)
  dxdzchi = dxdzchi_gs(il:iu, jl:ju, kl:ku)
  dydychi = dydychi_gs(il:iu, jl:ju, kl:ku)
  dydzchi = dydzchi_gs(il:iu, jl:ju, kl:ku)
  dzdzchi = dzdzchi_gs(il:iu, jl:ju, kl:ku)

  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              if (lapse) then
                 alp(i,j,k) = 1._dp
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
              
              ! g_ij = [1 - 2 R_c] \delta_ij + \partial_i \partial_j \chi
              gxx(i,j,k) = asq * (1._dp - 2._dp * rc(i,j,k) + dxdxchi(i,j,k))
              gxy(i,j,k) = asq * dxdychi(i,j,k)
              gxz(i,j,k) = asq * dxdzchi(i,j,k)
              gyy(i,j,k) = asq * (1._dp - 2._dp * rc(i,j,k) + dydychi(i,j,k))
              gyz(i,j,k) = asq * dydzchi(i,j,k)
              gzz(i,j,k) = asq * (1._dp - 2._dp * rc(i,j,k) + dzdzchi(i,j,k))

              ! K_ij = -a' [1 - 2 R_c] \delta_ij + ( a H'/H - a' ) \partial_i \partial_j \chi
              kxx(i,j,k) = -adot * gxx(i,j,k) + dxdxchi(i,j,k) * (a0 * hubdot / hub)
              kxy(i,j,k) = -adot * gxy(i,j,k) + dxdychi(i,j,k) * (a0 * hubdot / hub)
              kxz(i,j,k) = -adot * gxz(i,j,k) + dxdzchi(i,j,k) * (a0 * hubdot / hub)
              kyy(i,j,k) = -adot * gyy(i,j,k) + dydychi(i,j,k) * (a0 * hubdot / hub)
              kyz(i,j,k) = -adot * gyz(i,j,k) + dydzchi(i,j,k) * (a0 * hubdot / hub)
              kzz(i,j,k) = -adot * gzz(i,j,k) + dzdzchi(i,j,k) * (a0 * hubdot / hub)

              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! perturb the matter
                 press(i,j,k) = 0._dp ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0._dp
                 rho(i,j,k)   = rho0 * (1._dp + delta(i,j,k))
                 vel(i,j,k,:) = 0._dp
              endif

           endif

        enddo
     enddo
  enddo

  !
  ! make sure the metric type is physical (not conformal)
  !
  call check_metric()


end subroutine FLRW_SynchComoving
