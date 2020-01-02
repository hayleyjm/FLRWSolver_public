! file    flrw_powerspectrum.F90
! author  Hayley Macpherson
! date    30.12.2019
! desc    A spectrum of perturbations to FLRW initial data, longitudinal gauge, zero shift
!            read-in from files *only* at the moment -- IC generator to be incorporated! 

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
  CCTK_REAL :: a0,rho0,asq,rhostar,hub,adot,hubdot
  CCTK_REAL :: phi_ijk,alpijk

  ! globally-size arrays (to read in initial data files)
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3))   :: phi_gs, delta_gs
  CCTK_REAL, dimension(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3),3) :: delta_vel_gs

  ! locally-sized arrays (for this processor)
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))   :: phi, delta
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),3) :: delta_vel

  character(len=100) :: deltafile,vel1file,vel2file,vel3file,phifile,ics_dir
  integer :: dr_unit,dv_unit1,dv_unit2,dv_unit3,p_unit
  integer :: ics_dir_len,il,jl,kl,iu,ju,ku
  
  call CCTK_INFO("Initialising a POWER SPECTRUM of perturbations to an FLRW spacetime")
  
  !
  ! set logicals that tell us whether we want to use FLRWSolver to set ICs
  !
  call set_logicals(lapse,dtlapse,shift,data,hydro)

  !
  ! set parameters used in setting metric, matter parameters
  !
  call set_parameters(a0,rho0,asq,rhostar,hub,adot,hubdot)

  !
  ! read in perturbations from files
  !

  ! convert CCTK_STRING "FLRW_ICs_dir" to Fortran string
  call CCTK_FortranString(ics_dir_len,FLRW_ICs_dir,ics_dir)
 
  write(deltafile,'(a,a)')trim(ics_dir),'/delta.dat'
  write(vel1file,'(a,a)')trim(ics_dir),'/vel1.dat'
  write(vel2file,'(a,a)')trim(ics_dir),'/vel2.dat'
  write(vel3file,'(a,a)')trim(ics_dir),'/vel3.dat'
  write(phifile,'(a,a)')trim(ics_dir),'/phi.dat'

  open(newunit=dr_unit,file=deltafile,status='old')  ! delta rho file
  open(newunit=dv_unit1,file=vel1file,status='old')  ! delta vel file [1]
  open(newunit=dv_unit2,file=vel2file,status='old')  ! delta vel file [2]
  open(newunit=dv_unit3,file=vel3file,status='old')  ! delta vel file [3] 
  open(newunit=p_unit,file=phifile,status='old')     ! phi file

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
  delta     = delta_gs(il:iu, jl:ju, kl:ku)
  delta_vel = delta_vel_gs(il:iu, jl:ju, kl:ku, :)
  phi       = phi_gs(il:iu, jl:ju, kl:ku)

  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !
           ! set up metric, extrinsic curvature, lapse and shift
           !
           if (data) then

              phi_ijk = phi(i,j,k)
              
              if (lapse) then
                 alpijk     = FLRW_lapse_value * sqrt(1._dp + 2._dp * phi_ijk)
                 alp(i,j,k) = alpijk
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

              kxx(i,j,k) = -adot * a0 * (1._dp - 2._dp * phi_ijk) / alpijk
              kxy(i,j,k) = 0._dp
              kxz(i,j,k) = 0._dp
              kyy(i,j,k) = -adot * a0 * (1._dp - 2._dp * phi_ijk) / alpijk
              kyz(i,j,k) = 0._dp
              kzz(i,j,k) = -adot * a0 * (1._dp - 2._dp * phi_ijk) / alpijk

              !
              ! set up  matter variables
              !
              if (hydro) then
                 !
                 ! perturb the matter
                 press(i,j,k) = 0._dp ! pressure will be overwritten by EOS_Omni anyway
                 eps(i,j,k)   = 0._dp
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
