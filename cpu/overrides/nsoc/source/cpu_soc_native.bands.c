/* The unified target-band stage consumes H/S after the native force phase. */
#include <stdio.h>
#include "cpu_band_options.h"
#define CPU_BAND_USES_SPARSE_HS() (CpuBandOptions.band_on)
#include "../../soc/source/cpu_soc_native.tr.c"
