'''
#
# A script to take the annoyingly-formatted Cactus HDF5 3D data output and make it easier to use
#
#    Run me from your simulation directory.
#
#    Will take all HDF5 outputs in the directory (*.xyz.h5), and split these into a single file for each timestep
#         which contains all variables you have output. Nice.
#
#                Author: Hayley Macpherson (2018)
#
'''


import numpy as np
import h5py
import subprocess
import time
import os

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

# -----------------------------
# Name the new .hdf5 files based on the simulation name.
#
#    - get current directory name and delete everything before the LAST slash
#        using sed (this works for both dirs with _restart.par AND with simfactory-style directory structure)
# -----------------------------

# pwd = str(BASH('pwd'))                                    # Find directory we are in
# idx = pwd.find('FLRW')                                    # Find index where 'FLRW' starts in directory path
# dr = pwd[idx:]                                            # Extract directory from BASH command
#print 'Current directory is {}/'.format(dr)               # Confirm directory name

#par = str(BASH('ls *.par'))                                 # all simulation directories should include the parameter file with the same name (NOTE some will have _restart.par...)
#dr  = par[:-4]                                              # remove the .par at the end of the string to get sim name
#print(f" par = {par}, dr = {dr}")

pwd = str(BASH("pwd"))                   # Current directory
dr  = str(BASH("pwd | sed 's:.*/::'"))   # Strip the simulation name (same as the par file) from the current location

# check if any split files already exist - tell user to delete them as script will fail with files existing...
if os.path.isfile(f'./{dr}_it0.hdf5'):
    print('WARNING!!! Some split files may already exist: DELETE THEM BEFORE PROCEEDING!')

files = BASH('ls *.xyz.h5').split('\n')                   # Files as all HDF5 output by Cactus (excludes chkpts)
print('============================================')
print('')
print('             HDF5 FILE SPLITTER             ')
print('')
print('============================================')
print(f' splitting files for {dr}:')
for fl in files:
    print(fl)
print('')

itcnt = 1    # count the number of iterations
modit = 1    # only split every modit'th iteration
print('splitting every {} iterations'.format(modit))


fcnt = -1
for file in files:
    if 'Ricci' in file:
        print('Ignoring Ricci')
        continue            # ignore Ricci (from ADMAnalysis)
    fcnt+=1                                 # Counter for files
    f = h5py.File(file,'r') # open the file as an object
    print('--------------------------------------------')
    print(f' file {fcnt+1} of {len(files)} ')
    print(f' filename: {f.filename} ... ')
    start=time.time()
    try:
        filekeys = list(f.keys())
        h5ls = False
    except:
        # bash to get list of keys with h5ls
        filekeys = BASH('h5ls ' + file).split('\n')
        h5ls = True
    print(f"   got keys in {time.time()-start} sec ")
    
    it_prev = -1
    cnt = 0 # counter for number of keys per iteration
    
    for key in filekeys:
#    for key in f:
        keystr = str(key)

        if h5ls: # clean string up to look like output of f.keys()
            idx_endofkey = keystr.find('Data',0)                             # find where dataset=stuff starts
            keystr = keystr[:idx_endofkey-1]                                 # -1 to remove the space before dataset=stuff
#            keystr = str(BASH("echo '{}' | tr -d '\\'".format(keystr)))      # delete backslashes from string
#            keystr = str(BASH("echo '{}' | sed 's/\\//g'".format(keystr)))    # sed might remove warnings from tr
            keystr = str(BASH("echo '{}' | sed 's/[\]//g'".format(keystr)))

        it_start = keystr.find('it=',0) # first index of iteration string, if=-1, then continue loop past this key
        if (it_start==-1):
            continue
        else:
            it_end = keystr.find(' ',it_start+3) # last index of iteration, start looking at it_start for the next space
            it = keystr[it_start+3:it_end] # +2  on start for t=, then +1 to get to first number
            print(keystr)
            it = int(it)
            if (it%modit==0 or itcnt==1):
#                print key
                if (it==it_prev and cnt>0): # we're still within the same iteration, but not the first dataset
                    f.copy(keystr,fnew) # copy dataset to new file

                elif (it!=it_prev and cnt>0 and key!=filekeys[-1]): # we've reached a new iteration, set cnt=0 to open a new file!
                    fnew.close()
                    cnt = 0     # reset counter to open a new file
                    itcnt += 1  # count number of independent iterations we have

                if (cnt==0):    # we're at the first dataset for this iteration
                    fname = f'{dr}_it{it:06d}.hdf5'
                    fnew = h5py.File(fname,'a') # open a file and create first dataset (a appends to existing file 'add')
                
                    if (fcnt==0): # we are on the first file 
                        # extract parameters and grid structure from file, add in new "Datasets" with variable names

                        # do nothing for h5ls  - as params and global attrs doesn't come up in h5ls
                        if not h5ls:
                            filekeylast = filekeys[-1]
                            f.copy(filekeylast,fnew)

#                        grp = fnew.create_group(filekeys[-1]) # create group with same name as last key in file - params
#                        dset = grp.create_dataset("All Parameters",data=f.get(filekeys[-1]).get("All Parameters"))
#                        dset = grp.create_dataset("Grid Structure v5",data=f.get(filekeys[-1]).get("Grid Structure v5"))
#                        dset = grp.create_dataset("Datasets",data=find_var(files))
                
#                    f.copy(key,fnew) # copy dataset to new file
                    f.copy(keystr,fnew)

                cnt+=1               # iterate counter, counting number of dsets per iteration    
                it_prev = it         # set it_prev to be this iteration, before moving onto the next one!

    f.close() # close the file

print('Done!')
