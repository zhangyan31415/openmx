#ifndef CPU_VNA_STREAM_H
#define CPU_VNA_STREAM_H
#include <stddef.h>
/* Gradient has three Cartesian planes, each norb*projectors native FP32s.
 * The callback may weight/consume it before the next pair reuses the buffer. */
typedef void (*CpuVnaPairGradient)(void *context,int local_atom,int neighbor,
                                 float *gradient,int norb,int projectors);
double CpuVnaDeferred_Stream(float *****view,CpuVnaPairGradient consume,void *context);
/* Optional unified object; standalone SOC builds retain the native force path. */
int CpuVnaForceStream_Mode(void) __attribute__((weak));
int CpuVnaForceStream_Run(double *****density) __attribute__((weak));
#endif
