'''

 Script to read in CAMB power spectrum and make Gaussian random density, + corresponding velocity and potential (phi) fields
 Based on Equations (26) from Macpherson et al. 2017 (linear perturbations to FLRW)

    --> This version works embedded in FLRWSolver using CFFI (i.e. this module is "called" from Fortran via C...)

  Author: Hayley J. Macpherson
 Created: 07/12/2016
 Updated: 11/08/2021 (check repo for most recent change date)

'''
import time
import numpy as np
from scipy import fftpack
from scipy.interpolate import interp1d

def make_ics(a_init,Hini,box_size,bsize_code,resol,num_ghosts,rseed,ierrfile='ierrs.txt'):
    '''
    Make Gaussian random initial conditions for FLRWSolver (called from flrw_powerspectrum.F90)

    a_init     : initial FLRW background scale factor
    Hini       : initial FLRW background Hubble rate
    resol      : numerical resolution of the simulation
    num_ghosts : number of ghost cells in each dimension
    box_size   : physical size of each side of the box in comoving Mpc/h
    bsize_code : size of the box in code units
    rseed      : random seed to use to generate the Gaussian random field for density perturb

    ierrfile   : name of the text file to write error flags to that will then be checked with Fortran.
                 this is because CCTK doesn't exit nicely (and sometimes carries on) if something fails in this code

    '''
    #
    # Error flags are output into a file called create_ics.err and are for common failures:
    #    1. powerspectrum failed to load: pk_ierr = 0 (1) if ok (fail)
    #    2. writing ICs to files fails
    ifile = open(ierrfile,"a") # 'a' is append, since we already added a flag in builder.py

    #
    # Names of IC's files to be written and then read in by FLRWSolver
    #
    deltafile = "init_delta.dat"
    phifile   = "init_phi.dat"
    vel1file  = "init_vel1.dat"
    vel2file  = "init_vel2.dat"
    vel3file  = "init_vel3.dat"

    # A text file containing the name of the power spectrum file (incl. path)
    #    --> this file is MADE by FLRWSolver at initialisation, reading from the .par file...
    #
    #    * I hate this. But passing a C string from Fortran to Python proved to be a nightmare.
    #        and I couldn't quite figure out passing an array either. So, this will do for now...
    #
    pspecnamefile = "pk_filename.txt"

    # Read powerspectrum filename as text in pspecnamefile, then read data itself...
    pkfile               = open(pspecnamefile)
    matterpower_filename = pkfile.read().rstrip()
    pkfile.close()
    try:
        matterpower_file = np.loadtxt(matterpower_filename,skiprows=1)
        pk_ierr          = 0 # execution OK
    except:
        matterpower_file = np.zeros([10,2]) # dummy array
        pk_ierr          = 1 # execution FAILED
    ifile.write(f'{pk_ierr}\n')

    # Note k and P are dimensionfull
    k = matterpower_file[:,0]                 # |k| values (in h/Mpc)
    P = matterpower_file[:,1]                 # P(k) values (in (Mpc/h)^3)

    res     = resol + num_ghosts              # resolution including ghost cells
    box_res = [resol,resol,resol]             # dimensions of 3-D arrays

    bsize   = [box_size,box_size,box_size]    # dimension of box in Mpc/h
    spacing = bsize_code / resol              # grid spacing in code units
    spacing_phys = box_size / resol           # grid spacing in Mpc/h

    #
    # 0. Create kx,ky,kz arrays in physical units
    #     because we need modk in physical units to interpolate power spectrum
    kx1_phys    = 2. * np.pi * np.fft.fftfreq(resol,d=spacing_phys)
    kxp,kyp,kzp = np.meshgrid(kx1_phys,kx1_phys,kx1_phys)
    modk_phys   = np.sqrt(kxp**2 + kyp**2 + kzp**2)

    #
    # 1. Create kx,ky,kz arrays in code units for delta --> phi,vel
    kx1      = 2. * np.pi * np.fft.fftfreq(resol,d=spacing)
    kx,ky,kz = np.meshgrid(kx1,kx1,kx1)
    modk2    = kx**2 + ky**2 + kz**2                  # calculate modk^2 in the grid

    #
    # 2. Add in zero frequency to |k| and P(k)
    k_new = np.insert(k,0,0)                          # add zero frequency and power to arrays
    P_new = np.insert(P,0,0)
    # 3. Interpolate power spectrum
    P_interp = interp1d(k_new,P_new)                  # interpolate using multi-mode power spectrum

    #
    # 4. Calculate Gaussian random field (delta) from P(k)
    #
    # The below may fail if powerspectrum doesn't sample full k range needed -- could add a flag
    if (pk_ierr==0):
        # a) Generate random realisation in k-space
        Re_random_ft = np.random.normal(size=box_res) # real part
        Im_random_ft = np.random.normal(size=box_res) # imaginary part
        random_ft    = Re_random_ft + 1j*Im_random_ft

        # b) Define scale factor to make P(k) dimensionless
        volume = np.product(bsize)   # volume of the domain in Mpc/h
        dV     = spacing_phys**3     # volume element, i.e. dx^3, in Mpc/h
        scale  = dV**2/volume

        # Scale random to power spectrum
        Pkscale    = P_interp(modk_phys)/scale   # Scaled (dimensionless) power spectrum
        random_ft *= np.sqrt(Pkscale)            # Scale random field by sqrt(P)

        # Inverse FT and take real part to get random
        random = np.real(np.fft.ifftn(random_ft))
    else:
        # powerspectrum file wasn't read properly; set random to zero so we can carry on nicely
        random = np.zeros([resol,resol,resol])

    # Ensure delta is centered around zero (still follows P(k)) & transform to Fourier space
    delta_sync    = random - random.mean()
    delta_sync_ft = np.fft.fftn(delta_sync)

    # Define constants C1, C3 from Macpherson et al. 2017 in code units
    C1 = 2. / (3. * Hini**2)            # equiv to: a_init / ( 4. * np.pi * Grhostar)
    C3 = - 2. / (3. * a_init * Hini)    # equiv to: - np.sqrt( a_init / ( 6. * np.pi * Grhostar ) ) / a_init

    #
    # 5. a) calculate phi from synchronous density
    #        we have a k=0 point at [0,0,0] in modk2;
    #        this will give phi-->NaN,inf and give a warning when we divide by zero
    #        but for now we just set it to zero directly afterwards
    phi_ft           = np.zeros(box_res) + 1j * np.zeros(box_res)
    phi_ft[1:,1:,1:] = - delta_sync_ft[1:,1:,1:] / (C1 * modk2[1:,1:,1:])
    phi              = np.real(np.fft.ifftn(phi_ft))

    #
    # 5. Calculate delta_long & delta_vel from phi
    #
    delta_ft = - phi_ft * (C1 * modk2 + 2.)
    delta    = np.real(np.fft.ifftn(delta_ft))

    velx_ft = C3 * 1j * kx * phi_ft
    velx    = np.real(np.fft.ifftn(velx_ft))

    vely_ft = C3 * 1j * ky * phi_ft
    vely    = np.real(np.fft.ifftn(vely_ft))

    velz_ft = C3 * 1j * kz * phi_ft
    velz    = np.real(np.fft.ifftn(velz_ft))

    #
    # 5. Put data into arrays including ghost cells
    # (the ghost values are here set to zero. Cactus fills these out for us after FLRWSolver is called.)
    #
    delta_gh = np.zeros([res,res,res])
    phi_gh   = np.zeros([res,res,res])
    velx_gh  = np.zeros([res,res,res])
    vely_gh  = np.zeros([res,res,res])
    velz_gh  = np.zeros([res,res,res])

    delta_gh[3:-3,3:-3,3:-3] = delta
    phi_gh[3:-3,3:-3,3:-3]   = phi
    velx_gh[3:-3,3:-3,3:-3]  = velx
    vely_gh[3:-3,3:-3,3:-3]  = vely
    velz_gh[3:-3,3:-3,3:-3]  = velz

    #
    # Convert 3-D arrays to 2-D arrays (x-y slices at each z-value)
    #    and write data to files
    #
    delta_2d = np.zeros([res**2,res])
    phi_2d   = np.zeros([res**2,res])
    velx_2d  = np.zeros([res**2,res])           # arrays to hold data to write to files
    vely_2d  = np.zeros([res**2,res])
    velz_2d  = np.zeros([res**2,res])

    for i in range(0,res):
        delta_2d[i*res:(i+1)*res,:] = delta_gh[:,:,i]
        phi_2d[i*res:(i+1)*res,:]   = phi_gh[:,:,i]
        velx_2d[i*res:(i+1)*res,:]  = velx_gh[:,:,i]
        vely_2d[i*res:(i+1)*res,:]  = vely_gh[:,:,i]
        velz_2d[i*res:(i+1)*res,:]  = velz_gh[:,:,i]

    #
    # Open the files for writing ('w' overwrites files)
    #
    try:
        if (pk_ierr==0):
            # don't create files unless we have initial data
            delta_file = open(deltafile,"w")
            phi_file   = open(phifile,"w")
            vel1_file  = open(vel1file,"w")
            vel2_file  = open(vel2file,"w")
            vel3_file  = open(vel3file,"w")

            np.savetxt(delta_file,delta_2d)
            np.savetxt(phi_file,phi_2d)
            np.savetxt(vel1_file,velx_2d)
            np.savetxt(vel2_file,vely_2d)
            np.savetxt(vel3_file,velz_2d)

            delta_file.close()
            phi_file.close()
            vel1_file.close()
            vel2_file.close()
            vel3_file.close()
        file_ierr = 0 # execution OK
    except:
        file_ierr = 1 # execution FAILED
    ifile.write(f'{file_ierr}\n')
    ifile.close()
