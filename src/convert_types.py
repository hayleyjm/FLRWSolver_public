'''

   Desc: Convert C pointers to corresponding Python type
 Author: Hayley J. Macpherson
   Date: 10/08/2020

 Notes:
    --> from my_module.py in embedpython/testing/ directory, see this file for asarray and aschar (in progress)
    --> my_module was adapted from: https://www.noahbrenowitz.com/post/calling-fortran-from-python/

'''
import numpy as np

# Create the dictionary mapping ctypes to np dtypes.
ctype2dtype = {}

# Integer types
for prefix in ('int', 'uint'):
    for log_bytes in range(4):
        ctype = '%s%d_t' % (prefix, 8 * (2**log_bytes))
        dtype = '%s%d' % (prefix[0], 2**log_bytes)
        # print( ctype )
        # print( dtype )
        ctype2dtype[ctype] = np.dtype(dtype)

# Floating point types
ctype2dtype['float'] = np.dtype('f4')
ctype2dtype['double'] = np.dtype('f8')

#
# convert the pointer 'ptr' to its associated type in python
#

# a single number - real OR integer
def asnum(ffi, ptr, **kwargs):
    # Get the canonical C type of ptr as a string.
    T = ffi.getctype(ffi.typeof(ptr).item)
    
    if T not in ctype2dtype:
        raise RuntimeError("Cannot create an array for element type: %s" % T)

    # bit of a hack, but only for a single number and this is the easiest way I could figure out, i.e. [0] in array below
    a = np.frombuffer(ffi.buffer(ptr), ctype2dtype[T], count=1)[0]
    return a
