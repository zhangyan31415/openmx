#ifndef CPU_GSOC_COMPACT_H
#define CPU_GSOC_COMPACT_H
#include <stddef.h>
#include <stdint.h>
#include "mpi.h"
/* Include openmx_common.h first for dcomplex. Rows are zero-based spinor rows;
 * local logical bands are sorted by the same (energy, original index) key as
 * native GSOC. A TR partner has its own occupation but no second Q allocation. */
typedef struct {
 int n,ofs,dtype,partner,child,nr,nc,nb,prow,pcol,nprows,npcols;
 void *Q;
} CpuGSOCCompactSector;
int CpuGSOC_CompactEnabled(void);
int CpuGSOC_CompactLocalBands(void);
int CpuGSOC_CompactSortedBand(int local_band);
int CpuGSOC_CompactPack(int row0,int rows,int band0,int bands,void *out,int ld);
size_t CpuGSOC_CompactBytes(void);
void CpuGSOC_CompactRelease(void);
int CpuGSOC_CompactGetPhysical(int band1,dcomplex *out);
/* Build transfers ownership of independent input Q allocations, nulling each
 * pointer after its bounded redistribution. The caller keeps ownership of input
 * basis arrays; Build makes its own CSR copy, so caller may free them immediately.
 * The compact store survives CpuGSOC_Release; release only after BOTH DM and EDM.
 * All ranks call Build collectively. */
void CpuGSOC_CompactBuild(MPI_Comm comm,int ns,int nsec,
 CpuGSOCCompactSector *sec,const int64_t *cp,const int *ri,
 const dcomplex *cv,const double *unsorted_eval,int maxn,double *ko);
#endif
