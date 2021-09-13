#include "cctk.h"

#include <stdlib.h>

void CCTK_FNAME(FLRW_SetPythonModulePath)(void)
{
  setenv("FLRW_SRCDIR", FLRW_SRCDIR, 1);
}
