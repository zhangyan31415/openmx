#ifndef OPENMX_CPU_NSOC_SYM_H
#define OPENMX_CPU_NSOC_SYM_H
#include <stddef.h>
#include "cpu_nsoc_owner.h"
int CpuNSym_Enabled(void);
void CpuNSym_ReserveEigenScratch(double **coeff,size_t native_capacity,size_t matrix_capacity);
void CpuNSym_ReleaseEigenScratch(void);
size_t CpuNSym_EigenAllocationCapacity(size_t native_capacity);
int CpuNSym_NativeReturnEnabled(void);
void CpuNSym_SetNativeReturn(int step,double *coeff,size_t capacity,int first,int last);
int CpuNSym_NativeReturned(int step);
int CpuNSym_StateProvider(int n,int first,int last,CpuNORowProvider*provider);
int CpuNSym_StateCutNorms(int first,int last,int count,const int*ids,long double*norms);
int CpuNSym_StateProviderAudit(double*max_abs);
/* After all DM/EDM, tail and AfterDensity consumers: discard Q pages while
 * preserving the caller-owned allocation and fixed sparse-U row metadata. */
void CpuNSym_RetireStateCoefficients(int step);
int CpuNSym_FullProjected(void);
int CpuNSym_Prepare(int step);
int CpuNSym_RawSRebuild(void);
void CpuNSym_Solve(int step,int n,int maxn,double *ko,double *Ss,double *Hs,double *Cs,double *****nh,double ****olp);
int CpuNSym_ConvergenceBlocked(void);
void CpuNSym_RequestFull(void);
int CpuNSym_MixingAge(void);
int CpuNSym_FullSteps(void);
#endif
