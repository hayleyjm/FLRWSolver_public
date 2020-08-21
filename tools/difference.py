# SCRIPT TO FIND DIFFERENCE BETWEEN TWO FILES AND WRITE TO NEW FILE
# note: file length HAS to be the same. need some way to work around this

import numpy as np
import pylab as plt
import itertools

#sine_density = np.loadtxt("FLRW_sin_matchham/density_000000000.dat")
#ham_density = np.loadtxt("FLRW_ham10/density_000000000.dat")

file1 = raw_input("First file to compare: ")
print "Loading", file1
file_1 = np.loadtxt(file1)

file2 = raw_input("Second file to compare: ")
print "Loading", file2
file_2 = np.loadtxt(file2)

column1 = int(raw_input("Column of file 1  to compare: ")) - 1
column2 = int(raw_input("Column of file 2 to compare: ")) - 1
col_1 = file_1[:,column1]
col_2 = file_2[:,column2]

diff = col_1 - col_2
print"Difference is:",diff

x_axis = raw_input("Choose x or t: ")
if x_axis == 't':
    temp = np.arange(0,100001,2.4)
    x = 1 + np.sqrt(6*np.pi*1e-8)*temp
elif x_axis == 'x':
    x = np.arange(-240,241,24)

def writedat(filename, x, y, xprecision=5, yprecision=20):
    with open(filename,'w') as f:
        for a, b in itertools.izip(x, y):
            print >> f, "%.*g\t%.*g" % (xprecision, a, yprecision, b)

file_name = raw_input("Name of output file: ")
writedat(file_name,x,diff)
print "written to", file_name 
