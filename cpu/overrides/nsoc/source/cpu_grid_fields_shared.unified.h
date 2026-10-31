#ifndef CPU_GRID_FIELDS_SHARED_UNIFIED_H
#define CPU_GRID_FIELDS_SHARED_UNIFIED_H
#include <stddef.h>
#include <stdint.h>
#include <mpi.h>
/* One runtime layout shared by every translation unit. Scalar and spinor
 * runs retain their original packed 8- and 14-field payloads respectively. */
extern int cpu_grid_fields_width;
#define CPU_GF_DENSITY 0
#define CPU_GF_VXC (cpu_grid_fields_width)
#define CPU_GF_VPOT (2*cpu_grid_fields_width)
#define CPU_GF_REF (3*cpu_grid_fields_width)
#define CPU_GF_HARTREE (3*cpu_grid_fields_width+1)
#define CPU_GF_FIELDS (3*cpu_grid_fields_width+2)
#define CPU_GF_MAX_FIELDS 14
/* Collective, before allocation; rejects rank disagreement or live layout. */
int CpuGridFields_Configure(MPI_Comm world,int components);
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
