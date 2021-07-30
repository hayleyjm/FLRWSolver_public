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
extern void call_make_ics(double *a_init, double* rhostar, double *box_size, int32_t *resol, int32_t *num_ghosts, int32_t *rseed);
"""

#
# (import sys and update sys.path to include current dir seems necessary here
# maybe because of the way Python is initialised when this is called?
# just a raw import doesn't work)
#

# ==================================================================================
#        USER INFO:
#
# In "module", please change the THREE paths added to sys.path below to correspond to:
#
# 0, YOUR path to FLRWSolver src code (for, e.g., ``import create_ics'' statement)
# 1, YOUR path to c2raytools src code (for ``import c2raytools'' statement)
# 2, YOUR path to c2raytools MAIN src/ directory (for c2raytools __init__.py import statements)
#
# ==================================================================================

module = """
from ics_plugin import ffi
import sys
sys.path.insert(0, "/cosma/home/dp002/dc-macp1/software/flrwsolver/src/")
sys.path.insert(1, "/cosma/home/dp002/dc-macp1/software/flrwsolver/c2raytools3/src/")
sys.path.insert(2, "/cosma/home/dp002/dc-macp1/software/flrwsolver/c2raytools3/src/c2raytools3/")
import create_ics
import convert_types

@ffi.def_extern()
def call_make_ics(a_init,rhostar,box_size,resol,num_ghosts,rseed):

    a_init     = convert_types.asnum(ffi,a_init)
    rhostar    = convert_types.asnum(ffi,rhostar)
    box_size   = convert_types.asnum(ffi,box_size)
    resol      = convert_types.asnum(ffi,resol)
    num_ghosts = convert_types.asnum(ffi,num_ghosts)
    rseed      = convert_types.asnum(ffi,rseed)

    create_ics.make_ics(a_init,box_size,resol,num_ghosts,rseed)

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
