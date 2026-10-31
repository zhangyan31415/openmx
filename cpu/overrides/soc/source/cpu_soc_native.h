#ifndef CPU_SOC_NATIVE_H
#define CPU_SOC_NATIVE_H
int CpuSOCNative_Enabled(void);
void CpuSOCNative_RestartState(int ok);
void CpuSOCNative_MarkSolve(int step,double *****h,double *****im,double ****olp);
void CpuSOCNative_MarkSolved(int step);
void CpuSOCNative_BeforeExport(void);
#endif
