#ifndef CPU_NSOC_OWNER_OPENMX_H
#define CPU_NSOC_OWNER_OPENMX_H
int CpuNSOCOwnerMode(void);
void CpuNSOCOwnerValidate(int md_iter);
double CpuNSOCOwnerDensity(int n,int first,int last,int global_entries,const double*coeff,const double*eigenvalues,double*****dm,double*****edm,MPI_Comm world);
void CpuNSOCOwnerDestroy(void);
#endif
