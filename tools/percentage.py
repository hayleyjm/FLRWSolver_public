import numpy as np
import itertools

file1 = raw_input("Difference file: ")
print "Loading", file1
file_1 = np.loadtxt(file1)

file2 = raw_input("Numerical value file: ")
print "Loading", file2
file_2 = np.loadtxt(file2)

column1 = int(raw_input("Column of difference file: ")) - 1
column2 = int(raw_input("Column of numerical value file: ")) - 1
col_1 = file_1[:,column1]
col_2 = file_2[:,column2]

percentage = (col_1/col_2)*100

x_axis = raw_input("Choose x or t: ")
if x_axis == 't':
    x = np.arange(0,100001,2.4)
elif x_axis == 'x':
    x = np.arange(-240,241,24)

def writedat(filename, x, y, xprecision=5, yprecision=20):
    with open(filename,'w') as f:
        for a, b in itertools.izip(x, y):
            print >> f, "%.*g\t%.*g" % (xprecision, a, yprecision, b)

file_name = raw_input("Name of output file: ")
writedat(file_name,x,percentage)
print "written to", file_name
