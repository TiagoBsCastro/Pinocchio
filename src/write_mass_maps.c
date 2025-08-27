/* MASS_MAPS minimal diagnostics skeleton (clean rewrite) */
#include "pinocchio.h"
#ifdef MASS_MAPS
#include <math.h>
#include <stdint.h>

#ifndef PLC
#error "MASS_MAPS requires PLC"
#endif

#ifndef RECOMPUTE_DISPLACEMENTS
#error "MASS_MAPS requires RECOMPUTE_DISPLACEMENTS"
#endif

#ifndef HAVE_CFITSIO
#error "MASS_MAPS requires HAVE_CFITSIO"
#endif
#include <cfitsio.h>

#endif /* MASS_MAPS */
