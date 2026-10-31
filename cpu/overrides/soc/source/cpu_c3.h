#ifndef OPENMX_CPU_C3_H
#define OPENMX_CPU_C3_H
#include <stddef.h>
int CpuC3_Mode(void);
int CpuC3_WarmupActive(void);
void CpuC3_RequestFull(void);
int CpuC3_ConvergenceBlocked(void);
int CpuC3_FullSteps(void);
int CpuC3_MixingAge(void);
int CpuC3_CompactEnabled(void);
/* Local configuration; DFT validates MPI agreement before this call. */
void CpuC3_CompactBalanceConfigure(int sort,int cyclic);
int CpuC3_CompactLocalBands(void);
/* Local independent-sector coefficient payload, excluding small sparse plans. */
size_t CpuC3_CompactBytes(void);
int CpuC3_CompactSortedBand(int local);
int CpuC3_CompactPack(int row0,int rows,int band0,int bands,void *out,int ld);
void CpuC3_CompactRelease(void);
void CpuC3_ReleaseDenseOneShot(void);
void CpuC3_Solve(int step,int maxn,double *ko,double *****nh,double *****im,
                double ****overlap,int *is2,int *ie2,dcomplex *sf,
                dcomplex *h,dcomplex *q,dcomplex *ev);
#endif
