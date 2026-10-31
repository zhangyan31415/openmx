#ifndef CPU_GRID_FIELDS_SHARED_H
#define CPU_GRID_FIELDS_SHARED_H
#include <stddef.h>
#include <stdint.h>
#include <mpi.h>
enum { CPU_GF_DENSITY=0,CPU_GF_VXC=2,CPU_GF_VPOT=4,CPU_GF_REF=6,CPU_GF_HARTREE=7,CPU_GF_FIELDS=8 };
extern int cpu_grid_fields_active;
/* Must be called collectively at common truncation entry, before any branch.
 * Returns0/1 for agreed mode and -1 for invalid/mixed flags or enabled guard. */
int CpuGridFields_Agree(MPI_Comm world,int valid);
int CpuGridFields_Init(MPI_Comm world,int n1,int n2,int n3,size_t local_B,size_t local_C);
int CpuGridFields_Active(void);
int CpuGridFields_Requested(void);
double *CpuGridFields_Row(int field,int component);
int CpuGridFields_Owns(const double *ptr);
/* Return1 after a real owner-write publication,0 for an unregistered target,
 * -1 for a collective contract mismatch. These are never called inside omp. */
int CpuGridFields_Copy2(double **source,double **target,int components);
int CpuGridFields_Copy1(double *source,double *target);
void CpuGridFields_Destroy(void);
size_t CpuGridFields_Points(void);
size_t CpuGridFields_Bytes(void);
static inline size_t CpuGridFields_Index(size_t legacy_C,unsigned long long global_GN)
{return cpu_grid_fields_active?(size_t)global_GN:legacy_C;}
#endif
