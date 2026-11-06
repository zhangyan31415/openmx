/* Only the new H/S-only band path is independent of GSOC eigenvector buffers. */
#include <stdio.h>
#include "Inputtools.h"
#include "cpu_band_options.h"
#define CPU_BAND_USES_SPARSE_HS() (CpuBandOptions.band_on)
#include "../../soc/source/gsoc_runtime_guard.c"
