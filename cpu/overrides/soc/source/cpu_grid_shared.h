#ifndef OPENMX_CPU_GRID_SHARED_H
#define OPENMX_CPU_GRID_SHARED_H
#include <stddef.h>
/* Include mpi.h and define Type_Orbs_Grid (normally openmx_common.h) first. */
int CpuGridShared_Requested(void);
void CpuGridShared_Init(MPI_Comm comm, int natoms, const int *owners,
                        const int *grid_counts, const int *species,
                        const int *species_norb);
int CpuGridShared_IsAtom(int global_atom);
Type_Orbs_Grid *CpuGridShared_Row(int global_atom, int grid_index);
void CpuGridShared_Sync(void);
/* Canonical owner writes always use Row; ghost consumers bind ReadRow. */
Type_Orbs_Grid *CpuGridShared_ReadRow(int global_atom,int grid_index);
void CpuGridShared_Publish(void);
typedef struct {
  int enabled,local_numa,domains;
  size_t mapped_bytes,copy_bytes;
  unsigned long long epochs;
  double copy_rankmax_seconds,publish_rankmax_seconds;
} CpuGridSharedNumaStats;
void CpuGridShared_NumaStats(CpuGridSharedNumaStats *out);
void CpuGridShared_Destroy(void);
size_t CpuGridShared_Bytes(void);
/* Checked prefix-byte addition, also used by the allocation planner. */
int CpuGridShared_PlanBytes(size_t ngrid, size_t norb, size_t previous,
                            size_t *next);
#endif
