'''

A script to generate 3x different resolution initial conditions for the SAME distribution

    - based on a cut power spectrum as set in cut_powerspectrum.py
    - makes ICs for the lowest resolution based on this powerspectrum
    - then interpolates to get the two higher resolution ICs
    - make sure to change the sys.path.insert commands below to where your c2raytools3 is located

'''
import sys
import astropy.units as unit
import astropy.constants as const
from astropy.cosmology import WMAP9
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
#sys.path.insert(1, "/cosma/home/dp002/dc-macp1/software/flrwsolver/c2raytools3/src/")
#sys.path.insert(2, "/cosma/home/dp002/dc-macp1/software/flrwsolver/c2raytools3/src/c2raytools3/")
sys.path.insert(1, "/Users/hayleymac/Documents/codes/flrwsolver/c2raytools3/src/")
sys.path.insert(1, "/Users/hayleymac/Documents/codes/flrwsolver/c2raytools3/src/c2raytools3/")
import c2raytools3
from c2raytools3.power_spectrum import _get_dims, _get_k, power_spectrum_1d
import cut_powerspectrum
from cut_powerspectrum import boxL,pkcutfile,cutfac,res # res here is resolution we used to cut Pk

#
# Choose the resolutions you want
#   -- assumes [low,mid,high]
#
resolutions = [32,64,96]
#
# Number of ghost cells in each dimension (both sides combined)
num_ghosts  = 6

#
# Random seed for ICs generation
rseed = 4957398

#
# Initial scale factor
a_init = 1.0

#
# Dimensions of box in code units
#      (for interpolation -- this shouldn't actually change anything physical)
#
xmin = 0.
xmax = 1.

# -------------- end user changes ----------------------

#
# Get spacing, 1D x-values for each resolution
boxLcode = xmax - xmin

dx0 = boxLcode / resolutions[0]
dx1 = boxLcode / resolutions[1]
dx2 = boxLcode / resolutions[2]

x0 = np.arange(xmin,xmax,dx0)
x1 = np.arange(xmin,xmax,dx1)
x2 = np.arange(xmin,xmax,dx2)

print(f"       Box size = {boxL} cMpc")
#print(f"       Box size = {boxLcode} code units")
print(f"    Resolutions = {resolutions}")
print(f" Power spectrum = {pkcutfile}")


#
# Set some constants. We need length units in Mpc since our k's and P(k) are in units of Mpc^-1.
#
c = const.c.to('Mpc/s')                                  # speed of light in Mpc/s
G = const.G.to('Mpc^3/kg s^2')                           # gravitational constant in Mpc^3/(kg s^2)
rhostar = WMAP9.critical_density0.to('kg / Mpc^3')       # conserved FLRW density in kg/Mpc^3

C1 = a_init / ( 4. * np.pi * G * rhostar)                # constants C1, C3 from Macpherson et al. 2016
C3 = - np.sqrt( a_init / ( 6. * np.pi * G * rhostar ) )  # in physical units

icut = int(cutfac) # for file naming
#
# Names of IC's files (for each resolution}
#
deltafiles = []; phifiles = []
vel1files = []; vel2files = []; vel3files = []
for r in resolutions:
    deltafiles.append(f"init_delta_{r}_{icut}dx{res}.dat")
    phifiles.append(f"init_phi_{r}_{icut}dx{res}.dat")
    vel1files.append(f"init_vel1_{r}_{icut}dx{res}.dat")
    vel2files.append(f"init_vel2_{r}_{icut}dx{res}.dat")
    vel3files.append(f"init_vel3_{r}_{icut}dx{res}.dat")

#
# Make ICs for lowest resolution
#      (largely copied from create_ics.py)
#
try:
    matterpower_file = np.loadtxt(pkcutfile,skiprows=1)
except:
    print(f"ERROR opening file. Check your power spectrum: {pkcutfile}")

k = matterpower_file[:,0]          # |k| values
P = matterpower_file[:,1]          # P(k) values

resol   = resolutions[0]           # choose lowest res
box_res = [resol,resol,resol]      # dimensions of 3-D arrays

bsize   = [boxL,boxL,boxL]
spacing = boxL / resol             # grid spacing in cMpc

#
# Add in zero frequency to |k| and P(k)
#  & Interpolate power spectrum
#
k_new = np.insert(k,0,0)           # add zero frequency and power to arrays
P_new = np.insert(P,0,0)
P_interp = interp1d(k_new,P_new)   # interpolate using multi-mode power spectrum

#
# Make random field for lowest resolution
#
random = c2raytools3.gaussian_random_field.make_gaussian_random_field(box_res, boxL, P_interp, random_seed=rseed)
delta0  = random - random.mean() # Ensure delta is centered around zero (still follows P(k))

#
# Make an array which is size [res+1,res+1,res+1] so that extent is xmin,xmax INCLUSIVE
#    --> this means just implement periodic boundaries by addng repeat cell on the -1 side...
#
delta0_forinterp = np.zeros([resol+1,resol+1,resol+1])
delta0_forinterp[:-1,:-1,:-1] = delta0
delta0_forinterp[-1,:-1,:-1]  = delta0[0,:,:]
delta0_forinterp[:-1,-1,:-1]  = delta0[:,0,:]
delta0_forinterp[:-1,:-1,-1]  = delta0[:,:,0]
x0_forinterp = np.append(x0,x0[-1]+dx0) # add an extra point onto x-array to get xmax included
#
# Make an interpolating function
interp_delta0 = RegularGridInterpolator((x0_forinterp, x0_forinterp, x0_forinterp), delta0_forinterp)

#
# Interpolate lowest res to medium res
#
#      This is the form we need to pass x-points into RegGridInterp
pts1   = np.array(np.meshgrid(x1,x1,x1,indexing='ij')).T
delta1 = interp_delta0(pts1).T    # we need to Transpose this - from trial and error
delta1 = delta1 - delta1.mean()

#
# Set up for interpoation from medium to high res
#
resol = resolutions[1]
delta1_forinterp = np.zeros([resol+1,resol+1,resol+1])
delta1_forinterp[:-1,:-1,:-1] = delta1
delta1_forinterp[-1,:-1,:-1]  = delta1[0,:,:]
delta1_forinterp[:-1,-1,:-1]  = delta1[:,0,:]
delta1_forinterp[:-1,:-1,-1]  = delta1[:,:,0]
x1_forinterp = np.append(x1,x1[-1]+dx1)
#
# Make an interpolating function
interp_delta1 = RegularGridInterpolator((x1_forinterp, x1_forinterp, x1_forinterp), delta1_forinterp)

#
# Interpolate medium res to highest res
#
pts2   = np.array(np.meshgrid(x2,x2,x2,indexing='ij')).T
delta2 = interp_delta1(pts2).T
delta2 = delta2 - delta2.mean()

#
# Store delta's in a list to make looping easier
deltas = [delta0,delta1,delta2]
dxs    = [dx0,dx1,dx2]

#
# Loop over resolutions and make remaining ICs + write to files
#
for i,r in enumerate(resolutions):
    print(f" Making ICs for res = {r}")

    #
    # Make wavenumber arrays for FTs
    kx1      = 2. * np.pi * np.fft.fftfreq(r,d=dxs[i]) / unit.Mpc
    kx,ky,kz = np.meshgrid(kx1,kx1,kx1)
    modk2    = kx**2 + ky**2 + kz**2

    #
    # Move to FT and make phi,vel ICs
    print("    entering Fourier space ... ")
    delta_ft = np.fft.fftn(deltas[i])
    phi_ft   = - delta_ft / (C1 * modk2 + 2./c**2)
    phi      = np.real(np.fft.ifftn(phi_ft)) * phi_ft.unit

    velx_ft = C3 * 1j * kx * phi_ft
    velx    = np.real(np.fft.ifftn(velx_ft)) * velx_ft.unit

    vely_ft = C3 * 1j * ky * phi_ft
    vely    = np.real(np.fft.ifftn(vely_ft)) * vely_ft.unit

    velz_ft = C3 * 1j * kz * phi_ft
    velz    = np.real(np.fft.ifftn(velz_ft)) * velz_ft.unit
    print("     exiting Fourier space ... ")
    
    #
    # Convert to code units
    phi_on_c2 = phi / c**2
    velx_on_c = velx / c
    vely_on_c = vely / c
    velz_on_c = velz / c

    #
    # Store in arrays including ghost cells (set to zero)
    rgh      = r + num_ghosts
    delta_gh = np.zeros([rgh,rgh,rgh])
    phi_gh   = np.zeros([rgh,rgh,rgh])
    velx_gh  = np.zeros([rgh,rgh,rgh])
    vely_gh  = np.zeros([rgh,rgh,rgh])
    velz_gh  = np.zeros([rgh,rgh,rgh])

    delta_gh[3:-3,3:-3,3:-3] = deltas[i]
    phi_gh[3:-3,3:-3,3:-3]   = phi_on_c2
    velx_gh[3:-3,3:-3,3:-3]  = velx_on_c
    vely_gh[3:-3,3:-3,3:-3]  = vely_on_c
    velz_gh[3:-3,3:-3,3:-3]  = velz_on_c

    #
    # Write to 2D file structure we can read-in in FLRWSolver
    print(f"    writing to files ... ")
    delta_2d = np.zeros([rgh*rgh,rgh])
    phi_2d   = np.zeros([rgh*rgh,rgh])
    velx_2d  = np.zeros([rgh*rgh,rgh])
    vely_2d  = np.zeros([rgh*rgh,rgh])
    velz_2d  = np.zeros([rgh*rgh,rgh])
    for j in range(rgh):
        delta_2d[j*rgh:(j+1)*rgh,:] = delta_gh[:,:,j]
        phi_2d[j*rgh:(j+1)*rgh,:]   = phi_gh[:,:,j]
        velx_2d[j*rgh:(j+1)*rgh,:]  = velx_gh[:,:,j]
        vely_2d[j*rgh:(j+1)*rgh,:]  = vely_gh[:,:,j]
        velz_2d[j*rgh:(j+1)*rgh,:]  = velz_gh[:,:,j]
        
    delta_file = open(deltafiles[i],"w")
    phi_file   = open(phifiles[i],"w")
    vel1_file  = open(vel1files[i],"w")
    vel2_file  = open(vel2files[i],"w")
    vel3_file  = open(vel3files[i],"w")
    
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
print(" Done -- your initial conditions are in your current directory. ")
