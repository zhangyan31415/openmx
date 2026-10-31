#ifndef CPU_GSOC_H
#define CPU_GSOC_H
#include "cpu_gsoc_compact.h"
/* openmx_common.h must precede this header for dcomplex. */
int CpuGSOC_Enabled(void);
void CpuGSOC_RestartState(int ok);
void CpuGSOC_Solve(int step,int maxn,double*ko,double*****nh,double*****im,double****overlap,int*is2,int*ie2,dcomplex*sf,dcomplex*h,dcomplex*q,dcomplex*ev);
void CpuGSOC_Release(void);
#endif
