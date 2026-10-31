#ifndef CPU_SYMHISTORY_DIIS_H
#define CPU_SYMHISTORY_DIIS_H
#include <stddef.h>
#include <mpi.h>
typedef struct CpuSymHistoryDIIS CpuSymHistoryDIIS;
typedef void (*CpuSymHistoryInverse)(int n, double **a, double **ia, void *context);
typedef struct {
  int effective_step, num_mix, num_slide, fallback;
  double norm_current, norm_optimized, mixing_coefficient;
} CpuSymHistoryDIISInfo;
/* Owns 2*(history+2) compact vectors. No collective during creation. */
CpuSymHistoryDIIS *CpuSymHistoryDIIS_Create(size_t length, const double *weights,
                                         int history, MPI_Comm comm);
void CpuSymHistoryDIIS_Free(CpuSymHistoryDIIS *plan);
/* Local vector length may be zero. Slots 0..history+1. Returned pointers
 * must not be retained across Step(), which may rotate internal slots. */
double *CpuSymHistoryDIIS_Density(CpuSymHistoryDIIS *plan, int slot);
double *CpuSymHistoryDIIS_Residual(CpuSymHistoryDIIS *plan, int slot);
int CpuSymHistoryDIIS_Set(CpuSymHistoryDIIS *plan, int residual, int slot,
                         const double *values);
size_t CpuSymHistoryDIIS_Length(const CpuSymHistoryDIIS *plan);
/* First step may be native effective step2 or3; thereafter strictly sequential.
 * Adapter may call native branch2 and import D0,D1,D2,R2 before first step3.
 * Side effects on OpenMX NormRD/History_Uele are handled by native adapter.
 * Matrix inverse callback exactly follows native signature (n = augmented
 * order minus one) and should use the unchanged native Inverse algorithm.
 * Returns0 success,-1 invalid argument/order,-2 allocation,-3 nonfinite input.
 * All ranks must agree arguments and call order. Gram/Norm collectives use comm.
 */
int CpuSymHistoryDIIS_Step(CpuSymHistoryDIIS *plan, int effective_step,
                          double simple_weight, CpuSymHistoryInverse inverse,
                          void *context, CpuSymHistoryDIISInfo *info);
#endif
