'''

 Script to read in CAMB power spectrum and make Gaussian random density, + corresponding velocity and potential (phi) fields
 Based on Equations (26) from Macpherson et al. 2017 (linear perturbations to FLRW)

    --> This version works embedded in FLRWSolver using CFFI (i.e. this module is "called" from Fortran via C...)
    --> NOTE: changes to this file DO NOT require a Cactus re-compile!

  Author: Hayley J. Macpherson
 Created: 07/12/2016
 Updated: 10/08/2020 -- for compatability with CFFI embedding in FLRWSolver

'''
import time
import numpy as np
import c2raytools3 # there is a local (.tar) copy of this included in the repo
import astropy.units as unit
import astropy.constants as const
from scipy import fftpack
from c2raytools3.power_spectrum import _get_dims, _get_k, power_spectrum_1d
from scipy.interpolate import interp1d
from astropy.cosmology import WMAP9

def make_ics(a_init,box_size,resol,num_ghosts,rseed):
    '''
    Make Gaussian random initial conditions for FLRWSolver (called from FLRWSolver)

    a_init     : initial FLRW background scale factor
    resol      : numerical resolution of the simulation
    num_ghosts : number of ghost cells in each dimension
    box_size   : physical size of each side of the box in cMpc
    rseed      : random seed to use to generate the Gaussian random field for density perturb

    '''
    
    #
    # Set some constants. We need length units in Mpc since our k's and P(k) are in units of Mpc^-1.
    #
    c = const.c.to('Mpc/s')                                  # speed of light in Mpc/s
    G = const.G.to('Mpc^3/kg s^2')                           # gravitational constant in Mpc^3/(kg s^2)
    rhostar = WMAP9.critical_density0.to('kg / Mpc^3')       # conserved FLRW density in kg/Mpc^3

    C1 = a_init / ( 4. * np.pi * G * rhostar)                # constants C1, C3 from Macpherson et al. 2016
    C3 = - np.sqrt( a_init / ( 6. * np.pi * G * rhostar ) )  # in physical units

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
    except:
        print(f"ERROR: opening file. Check your power spectrum is really at {matterpower_filename}")

    k = matterpower_file[:,0]                 # |k| values
    P = matterpower_file[:,1]                 # P(k) values

    res     = resol + num_ghosts              # resolution including ghost cells
    box_res = [resol,resol,resol]             # dimensions of 3-D arrays

    bsize   = [box_size,box_size,box_size]
    spacing = box_size / resol                # grid spacing in cMpc

    #
    # 1. Create kx,ky,kz arrays
    # 2. Add in zero frequency to |k| and P(k)
    # 3. Interpolate power spectrum
    #
    kx1      = 2. * np.pi * np.fft.fftfreq(resol,d=spacing) / unit.Mpc
    kx,ky,kz = np.meshgrid(kx1,kx1,kx1)
    modk2    = kx**2 + ky**2 + kz**2                  # calculate modk^2 in the grid

    k_new = np.insert(k,0,0)                          # add zero frequency and power to arrays
    P_new = np.insert(P,0,0)
    P_interp = interp1d(k_new,P_new)                  # interpolate using multi-mode power spectrum

    #
    # 1. Calculate Gaussian random field (delta) from P(k)
    # 2. Calculate phi & delta_vel using FFT
    # 3. Include zero value ghost cells
    #

    # Arrays including ghost cells
    delta_gh = np.zeros([res,res,res])
    phi_gh   = np.zeros([res,res,res])
    velx_gh  = np.zeros([res,res,res])
    vely_gh  = np.zeros([res,res,res])
    velz_gh  = np.zeros([res,res,res])

    random = c2raytools3.gaussian_random_field.make_gaussian_random_field(box_res, box_size, P_interp, random_seed=rseed)
    # Ensure delta is centered around zero (still follows P(k))
    delta  = random - random.mean()

    delta_ft = np.fft.fftn(delta)
    phi_ft   = - delta_ft / (C1 * modk2 + 2./c**2)
    phi      = np.real(np.fft.ifftn(phi_ft)) * phi_ft.unit

    velx_ft = C3 * 1j * kx * phi_ft
    velx    = np.real(np.fft.ifftn(velx_ft)) * velx_ft.unit

    vely_ft = C3 * 1j * ky * phi_ft
    vely    = np.real(np.fft.ifftn(vely_ft)) * vely_ft.unit

    velz_ft = C3 * 1j * kz * phi_ft
    velz    = np.real(np.fft.ifftn(velz_ft)) * velz_ft.unit

    # Convert to code units (delta already dimensionless)
    phi_on_c2 = phi / c**2
    velx_on_c = velx / c
    vely_on_c = vely / c
    velz_on_c = velz / c

    #
    # Put data into arrays including ghost cells
    # (the ghost values are here set to zero. Cactus fills these out for us after FLRWSolver.)
    #
    delta_gh[3:-3,3:-3,3:-3] = random
    phi_gh[3:-3,3:-3,3:-3]   = phi_on_c2
    velx_gh[3:-3,3:-3,3:-3]  = velx_on_c
    vely_gh[3:-3,3:-3,3:-3]  = vely_on_c
    velz_gh[3:-3,3:-3,3:-3]  = velz_on_c

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

    np.savetxt(deltafile,delta_2d)
    np.savetxt(phifile,phi_2d)
    np.savetxt(vel1file,velx_2d)
    np.savetxt(vel2file,vely_2d)
    np.savetxt(vel3file,velz_2d)
