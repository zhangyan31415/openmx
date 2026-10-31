#ifndef CPU_FORCE_DM_SHARED_H
#define CPU_FORCE_DM_SHARED_H
#include <mpi.h>
#include <stddef.h>
typedef struct CpuForceDM CpuForceDM;
/* Inputs use native atom (1-based) and neighbor (0-based) order. All six
 * source components are local owned rows; no source ghost reads occur. */
CpuForceDM *CpuForceDM_Create(MPI_Comm comm, int atoms, int owned, int ghosts,
 const int *owner, const int *m2g, const int *fm2g, const int *norb,
 const int *fnan, int *const *natn, double *****cdm, double *****idm,
 double ******cdm_view, double ******idm_view);
void CpuForceDM_Destroy(CpuForceDM **cache);
size_t CpuForceDM_Bytes(const CpuForceDM *cache);
/* Pure checked byte planning, exposed for boundary tests. */
int CpuForceDM_AddBytes(size_t rows,size_t cols,size_t previous,size_t *next);
int CpuForceDM_Requested(MPI_Comm comm);
#endif
