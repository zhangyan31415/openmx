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
/* Collective fixed-geometry lifetime handoff. Suspend punches holes only in
 * complete payload pages of writable shared mappings (MADV_REMOVE), retaining
 * every window, offset, facade pointer and virtual capacity. Return 1 if any
 * window was discarded, 0 for a safe no-op. Partial platform support is legal:
 * regardless of the result all owners must call BeginRebuild, rewrite their
 * complete native orbital payload and Publish before any ReadRow consumer.
 * Direct Row is the owner-write API and remains address-stable throughout. */
int CpuGridShared_Suspend(void);
void CpuGridShared_BeginRebuild(void);
int CpuGridShared_IsSuspended(void);
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
