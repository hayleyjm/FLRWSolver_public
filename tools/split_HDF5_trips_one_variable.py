'''
#
# A script to take the mescaline trip reports and split into files containing only one chosen variable
#     -- this is useful if you want to share data with people but they only really care about one variable...
#     -- probably only for one or two circumstances
#
#    Run me from your simulation directory.
#
#                Author: Hayley Macpherson (2021)
#
'''


import numpy as np
import h5py
import subprocess
import time
import os

### define the variable you want
varsplit = 'rho_Aij'

##########################################################################################################
'''
Function to display keys for a file. NOTE: f.keys() *does* work, but this is a workaround for Python 3... don't ask
'''

def keys(f):
    return [key for key in list(f.keys())]


'''
Function to execute bash commands
*** BE CAREFUL RUNNING THIS AS COMMANDS INPUT ARE *ACTUALLY* EXCECUTED... SO DON'T DELETE ANYTHING!**
'''
def BASH(command):
    return subprocess.check_output(command,shell=True).decode().strip()


'''
Function to extract the name of the variable contained in a filename (Cactus HDF5 output)
	file :: can be a single filename, or a list of filenames
	returns :: list of, or single strings containing variable names

'''
def find_var(file):
    if (type(file)==list): # we have more than one file
        varnames = [] # empty list of var names
        for fl in file:
            f = h5py.File(fl,'r') # open file as an object
            fkeys0 = str(keys(f)[0])
            it_end = fkeys0.find(' ',0) # find the first space: this is where the variable name ends
            varnames.append(fkeys0[0:it_end])
            f.close()
        return varnames
            
    else: # we have one file
        f = h5py.File(file,'r')
        fkeys0 = str(keys(f)[0])
        it_end = fkeys0.find(' ',0) # find the first space: this is where the variable name ends
        f.close()
        return fkeys0[0:it_end]

##########################################################################################################

#
# the base of the files we are using
#     -- will use all files with this name in current directory
fbase = 'trip_report_3D'

#
# the name of the NEW files that will be created for each iteration
fbasenew = fbase + '_' + varsplit

files = BASH('ls ' + fbase + '*').split('\n')                   # Files as all HDF5 output by Cactus (excludes chkpts)
print('============================================')
print('')
print('             HDF5 FILE SPLITTER             ')
print('         (single-variable version)          ')
print('============================================')
print(f' splitting {varsplit} out of {fbase}*:')
for fl in files:
    print(fl)
print('')

itcnt = 1    # count the number of iterations

fcnt = -1
for file in files:
    fcnt+=1                                 # Counter for files
    f = h5py.File(file,'r') # open the file as an object
    print('--------------------------------------------')
    print(f' file {fcnt+1} of {len(files)} ')
    print(f' filename: {f.filename} ... ')
    
    # make the new filename with fbasenew for copying our variable to
    fnamenew = f.filename.replace(fbase,fbasenew)
    fnew     = h5py.File(fnamenew,'a') # open a file and create first dataset (a appends to existing file 'add')
    print(f' NEW file: {fnew.filename} ... ')
    
    # loop over keys and find the one we want, copy it to fnew
    filekeys = list(f.keys())
    for key in filekeys:
        keystr = str(key)
        # check if we're at the keyname we want
        if (key==varsplit):
            print(f'Got {key}, copying...')
            f.copy(key,fnew)

#                cnt+=1               # iterate counter, counting number of dsets per iteration
#                it_prev = it         # set it_prev to be this iteration, before moving onto the next one!

    f.close() # close the file
    fnew.close()
print('Done!')
