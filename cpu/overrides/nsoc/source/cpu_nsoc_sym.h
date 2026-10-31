#ifndef OPENMX_CPU_NSOC_SYM_H
#define OPENMX_CPU_NSOC_SYM_H
#include <stddef.h>
int CpuNSym_Enabled(void);
void CpuNSym_ReserveEigenScratch(double **coeff,size_t native_capacity,size_t matrix_capacity);
void CpuNSym_ReleaseEigenScratch(void);
int CpuNSym_FullProjected(void);
int CpuNSym_Prepare(int step);
int CpuNSym_RawSRebuild(void);
void CpuNSym_Solve(int step,int n,int maxn,double *ko,double *Ss,double *Hs,double *Cs,double *****nh,double ****olp);
int CpuNSym_ConvergenceBlocked(void);
void CpuNSym_RequestFull(void);
int CpuNSym_MixingAge(void);
int CpuNSym_FullSteps(void);
#endif
