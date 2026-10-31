#ifndef OPENMX_NATIVE_PANEL_SOLVER_H
#define OPENMX_NATIVE_PANEL_SOLVER_H
#include <stddef.h>
void CpuNSym_DirectConfigure(void);
int CpuNSym_DirectRegister(int step,int n,double *****h,double ****s);
void CpuNSym_DirectWorkspaceRecord(double *work,size_t nw,double *cdm,size_t nc,
 double *edm,size_t ne,double *pdm,size_t np);
void CpuNSym_DirectWorkspaceProve(double *work,double *cdm,double *edm,double *pdm);
#endif
