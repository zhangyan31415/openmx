#ifndef OPENMX_CPU_ORBITAL_WORKSPACE_H
#define OPENMX_CPU_ORBITAL_WORKSPACE_H
#include <stddef.h>
typedef struct CpuOrbitalWorkspace {
 int maxl,radial,active;size_t bytes,allocations;
 double **RF,**dRF,**AF,**dAFQ,**dAFP,**rows,*values;
} CpuOrbitalWorkspace;
CpuOrbitalWorkspace *CpuOrbitalWorkspace_Acquire(int maxl,int radial,double x,double y,double z,double cutoff);
void CpuOrbitalWorkspace_Release(CpuOrbitalWorkspace *workspace);
void CpuOrbitalWorkspace_ReleaseThread(void);
size_t CpuOrbitalWorkspace_ThreadAllocationCount(void);
#endif
