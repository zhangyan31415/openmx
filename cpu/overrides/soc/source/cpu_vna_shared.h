#ifndef CPU_VNA_SHARED_H
#define CPU_VNA_SHARED_H
#include <stddef.h>
#include <mpi.h>
/* Per-physical-node shared cache: own atoms plus actual FNAN atom readers.
 * FP32 ABI matches native Type_DS_VNA. Metadata is copied at initialization;
 * caller-owned atom/neighbor/image arrays may later be freed independently. */
int CpuVnaShared_Requested(void);
/* Call collectively and unconditionally before any flag-dependent Init. */
int CpuVnaShared_Agree(MPI_Comm comm);
int CpuVnaShared_Active(void);
void CpuVnaShared_ConfigureDeferred(int enabled);
int CpuVnaShared_Deferred(void);
/* Collective lifetime transition for the deferred, fixed-geometry path.
 * Drop invalidates raw values until owners rewrite them and Restore republishes
 * them to every needed node. Returned bytes are successful local-node-leader
 * MADV_REMOVE ranges, not a measured physical-memory saving. Unsupported shared
 * mappings retain their pages and report that fact; numerical behavior is safe. */
size_t CpuVnaShared_DropRaw0(void);
int CpuVnaShared_Raw0Dropped(void);
void CpuVnaShared_RestoreRaw0(void);
int CpuVnaShared_WorthDeferring(size_t local_dead_dense_bytes);
void CpuVnaShared_DerivativesReady(void);
double CpuVnaDeferred_Complete(float *****view);
void CpuVnaShared_Init(MPI_Comm comm,int atoms,const int*owners,const int*fan,
 int**natn,int**ncn,int**image4,const int*norb,int projectors,int scratch_fan,int scratch_norb);
float*CpuVnaShared_Row(int component,int atom,int neighbor,int orbital);
/* Collective world publication to all needed nodes.
 * Epoch0: raw values and derivatives; epoch1: derivatives weighted by ene. */
void CpuVnaShared_Publish(int epoch);
/* Collective world barrier: all raw readers complete before owners modify
 * derivative components. */
void CpuVnaShared_BeginWeighted(void);
void CpuVnaShared_BindAtom(float*****view,int scratch,int atom);
void CpuVnaShared_BindCenter(float*****view,int scratch,int center,int**reverse);
/* Actual shared allocation contributed by this rank (node leader only). */
size_t CpuVnaShared_OwnedBytes(void);
void CpuVnaShared_Destroy(void);
#endif
