import numpy as np
import pylab as plt

file1 = np.loadtxt('FLRW_alpha/hydrobase::rho.maximum.asc')
file2 = np.loadtxt('FLRW_alpha/scale_FLRW.out')

file1_col = file1[:,2]
file2_col = file2[:,3]
time = file1[:,1]

plt.close()
plt.plot(time,file1_col)
plt.plot(time,file2_col)
plt.xscale('log')
plt.xlabel('log time')
plt.ylabel('density')
plt.title('FLRW density evolution')
plt.show()
