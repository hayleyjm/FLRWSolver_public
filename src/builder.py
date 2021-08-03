'''

    Desc: Builds the static/dynamic library that links the Python create_ics code to FLRWSolver via C using CFFI
  Author: Hayley J. Macpherson (adapted from: https://www.noahbrenowitz.com/post/calling-fortran-from-python/)
    Date: 10/08/2020

  Details:
    --> this code makes a .dylib (Mac) or .so (Linux) file that can then be dynamically linked when compiling
         (note: I couldn't make this option work with the ET, you will need to figure it out)
    --> OR a C code that can be statically linked

'''

import cffi
ffibuilder = cffi.FFI()

static = True # if False; make dynamic library
if (static==False):
    print("WARNING: the current version of FLRWSolver is not written for dynamic library linking.")
    print("   either figure it out yourself, or set static = True.")

#
# NOTE the * in definitions here means we are defining a POINTER - which is what we want
# these particular integer types correspond to the DEFAULT Cactus integer type.
#
header = """
extern void call_make_ics(double *a_init, double *rhostar, double *box_size, double *bsize_code, int32_t *resol, int32_t *num_ghosts, int32_t *rseed);
"""

#
# (import sys and update sys.path to include current dir seems necessary here
# maybe because of the way Python is initialised when this is called?
# just a raw import doesn't work)
#

# ==================================================================================
#        USER INFO:
#
# In "module", please change 'flrwsolverpath' to correspond to where you have installed the FLRWSolver top directory
#
# Note:
#    1. Ensure you INCLUDE the trailing slash '/' in the FLRWSolver path
#    2. We assume you have used the un-tarred c2raytools3 included in this repo
#
# ==================================================================================

module = """
from ics_plugin import ffi
import sys
# --------- change path below ------------
flrwsolverpath="/path/to/flrwsolver/"
# ----------------------------------------

sys.path.insert(0, flrwsolverpath + "src/")
sys.path.insert(1, flrwsolverpath + "c2raytools3/src/")
sys.path.insert(2, flrwsolverpath + "c2raytools3/src/c2raytools3/")
import create_ics
import convert_types

@ffi.def_extern()
def call_make_ics(a_init,rhostar,box_size,bsize_code,resol,num_ghosts,rseed):

    a_init     = convert_types.asnum(ffi,a_init)
    rhostar    = convert_types.asnum(ffi,rhostar)
    box_size   = convert_types.asnum(ffi,box_size)
    bsize_code = convert_types.asnum(ffi,bsize_code)
    resol      = convert_types.asnum(ffi,resol)
    num_ghosts = convert_types.asnum(ffi,num_ghosts)
    rseed      = convert_types.asnum(ffi,rseed)

    create_ics.make_ics(a_init,rhostar,box_size,bsize_code,resol,num_ghosts,rseed)

"""

with open("plugin.h", "w") as f:
    f.write(header)

ffibuilder.embedding_api(header)
ffibuilder.set_source("ics_plugin", r'''
    #include "plugin.h"
''')

ffibuilder.embedding_init_code(module)
if (static):
    # Output the C code to be statically linked
    ffibuilder.emit_c_code("pspec_ics.c")
else:
    # Output a dynamic library to be dynamically linked...
    ffibuilder.compile(target="libplugin.*", verbose=True)
