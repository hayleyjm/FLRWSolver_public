'''

A script to cut a matter powerspectrum file output from CAMB
(or one in the same format)
to the desired minimum wavelength (maximum k)
        -- note most of this just cut + paste from old create_ics.py

Author: Hayley Macpherson
  Date: 24/08/2021

'''
import numpy as np

# ------------------------------------------------------
#
# User to change these settings below
#

# Physical box size (Mpc/h), initial redshift and resolution (number of grid cells in each dimension)
boxL = 1000.0
zini = 1000.0
res  = 32

# Name + location of power spectrum file
pspecfile = f"/Users/hayleymac/Documents/codes/flrwsolver/powerspectra/camb/FLRW_matterpower_z{int(zini)}.dat"
# Number of rows to skip, if you have a header in the file set this to 1
nskip = 1

# Factor to cut powerspectrum by as multiple of grid spacing
#      --> minimum wavelength = cutfac * dx = 2pi/kmax
cutfac = 10.0

# ------------------ end user changes ------------------

print(f"  ~~~ CUTTING POWER SPECTRUM ~~~   ")
print("")
print(f"    box length = {boxL} cMpc ")
print(f"    resolution = {res} ")
print("")
print("   ~~~                       ~~~ ")

# Check cutfac vs. Nyquist freq
if (cutfac==2):
    print(f"WARNING: setting cutfac = {cutfac} will do nothing (since this is Nyquist freq).")
elif (cutfac<2):
    print(f"WARNING: cannot set cutfac = {cutfac}, since this is below Nyquist freq.")
else:
    print(f"Cutting powerspectrum to a minimum wavelength of {cutfac} x dx")

# Load Pk data + cut
print(f"Loading powerspectrum data from {pspecfile}")
try:
    pdata = np.loadtxt(pspecfile,skiprows=nskip)
except:
    print(f"Error loading file at: {pspecfile}")
k = pdata[:,0]                 # |k| values
P = pdata[:,1]                 # P(k) values

P_cut  = np.zeros(P.size)                 # array for "cut" power spectrum
wavmin = cutfac * (boxL / res)            # minimum wavelength
kmax   = np.sqrt(3) * 2. * np.pi / wavmin # maximum k-value -- sqrt(3) from assuming kx=ky=kz

tol = 5e-02                               # tolerance for finding k-value in k array
try:
    index = np.where(np.isclose(k,kmax,rtol=tol)==True)[0][0]
except:
    # above will fail if tol is too small. increase it!
    print(f"Cutting P(k) failed with a tol = {tol}. Using tol = {10.*tol} instead.")
    index = np.where(np.isclose(k,kmax,rtol=10.*tol)==True)[0][0]
P_cut[:index+1] = P[:index+1] # +1 to make sure we include this value

pkcutfile = pspecfile + ".cut"
cutfile   = open(pkcutfile,"w")
cutdata   = np.zeros(pdata.shape)
cutdata[:,0] = k
cutdata[:,1] = P_cut
np.savetxt(cutfile,cutdata,header=f"k/h, P(k; min wavelength {cutfac} x dx)")
cutfile.close()
print("Done.")
print(f"Find your cut powerspectrum at {pkcutfile}. BYE")
