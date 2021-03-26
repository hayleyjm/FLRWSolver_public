#!/usr/bin/python
'''

A script to get the proper initial HL for use with FLRWSolver
    -- to make sure things are consistent with e.g. the box length
    
    USAGE:
    python get_init_HL.py res boxL_z0(Gpc) zini

    (mainly adapted from my scribblings in Projects/ImprovedAveraging/translate_phys_units.ipynb)
'''
import numpy as np
from astropy import units
from astropy import constants as const
import sys

'''
Define some functions for FLRW-stuff
'''
def xi(eta,rhostar,ainit,etainit):
    return 1 + np.sqrt(2.*np.pi*rhostar/(3.*ainit)) * (eta - etainit)

def aflrw(xival,ainit):
    return ainit * xival**2

def rhoflrw(xival,rhostar,ainit):
    return rhostar / aflrw(xival,ainit)**3

def hubflrw(eta,rhostar,ainit,etainit):
    xival = xi(eta,rhostar,ainit,etainit)
    return aflrw(xival,ainit) * np.sqrt(8.*np.pi*rhoflrw(xival,rhostar,ainit)/3.)

def get_z(aval,zinit):
    return (1.+zinit)/aval - 1.
    
    
'''
Define the initial params you want for the simulation
    -- use sys.argv, note sys.argv[0] is the name of script
====================================================================================
'''
# Proper/comoving length of box at redshift z=0
Lz0 = float(sys.argv[2]) * units.Gpc
#
# Initial redshift, scale factor, time (latter should stay the same)
zini  = float(sys.argv[3])
ainit = 1.0
# tinit = 1.0
#
# Simulation (code units) box size, dtfac, res, etc
res   = int(sys.argv[1])
boxL  = 1.0
dtfac = 0.1
dx    = boxL / res
dt    = dtfac * dx
#
# Redshift after which you'd like to increase freq of 3D output
zinc = 3.0
'''
====================================================================================
'''

#
# User info
print(f' Hello! finding initial HL for: ')
print('')
print(f'                   res = {res}')
print(f'    box length (z = 0) = {Lz0}')
print(f'    box length  (code) = {boxL}')
print(f'       initial reshift = {zini}')

print(f'')

'''
Assume EdS model -- i.e. H_0 ~ 45 km/s Mpc
    & find HL at z=0
'''
H0    = 45. * units.km / (units.s * units.Mpc)
dH_z0 = (const.c/H0).to('Mpc')
HLz0  = (H0*Lz0/const.c).to('')

print(f'                   H_0 = {H0}')
print(f'  Hubble horizon (z=0) = {dH_z0}')
print(f'               HL(z=0) = {HLz0}')
'''
Scale HL(z=0) back to desired initial redshift
'''
HL_zini = HLz0 * np.sqrt(1. + zini)
Lz0mpc = Lz0.to('Mpc')

'''
Find settings for final time, etc
    -- assume we want to run to z=0
    -- set up array of conformal times, translate to a_flrw
    -- then we can find final_time in conformal time to run to
'''
afinal      = ainit + zini # final scale factor we want to run to
Hinit       = HL_zini / boxL
tinit       = 2./(3.*Hinit)
print(f'      tinit = 2/3Hinit = {tinit}')
rhostarinit = Hinit**2 * 3. * ainit / (8.*np.pi)
etatest     = np.arange(tinit,1e5,dt)
xitest      = xi(etatest,rhostarinit,ainit,tinit)
aflrwval    = aflrw(xitest,ainit)
zvals       = get_z(aflrwval,zini)
aidx_fin    = np.where(aflrwval>afinal)[0][0]
zidx_inc    = np.where(zvals<zinc)[0][0] # index where z<1 for the first time, when to increase output
eta_inc     = etatest[zidx_inc]          # time at which we want to then increase freq. of 3D data
etafinal    = etatest[aidx_fin]
itfinal     = (etafinal-tinit)/dt
itz1        = (etatest[zidx_inc]-tinit)/dt  # iteration where z=1
print(f' running from a = {ainit} to a = {afinal} will take {int(itfinal)} iterations')
print(f'    and FYI you will reach z={zvals[zidx_inc]:.4f} at eta = {etatest[zidx_inc]:.4f} after {int(itz1)} iterations')


print(f'')
print(f'    ---> Settings for par file: ')
print(f' FLRWSolver::FLRW_init_HL            = {HL_zini}')
print(f' FLRWSolver::FLRW_init_a             = {ainit}')
print(f' FLRWSolver::FLRW_lapse_value        = {ainit}')
print(f' FLRWSolver::FLRW_boxlength          = {Lz0mpc}')
if (Lz0.value>2.0):
    print("WARNING: Make sure to use powerspectrum in LONGITUDINAL gauge for larger-scale sims")
else:
    print(f' FLRWSolver::FLRW_powerspectrum_file = FLRW_matterpower_z{int(zini)}.dat')

print('')
print(f'First run to z = {zinc}:')
print('')
print(f' Cactus::cctk_initial_time = {tinit}')
print(f' Cactus::cctk_final_time   = {eta_inc}')
print(f'   time::dtfac             = {dtfac} ')

print('')
print('--- RESTART PARAMS ---')
print(f' Cactus::cctk_final_time   = {etafinal}')
