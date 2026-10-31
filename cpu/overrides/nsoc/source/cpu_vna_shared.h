#ifndef CPU_VNA_SHARED_H
#define CPU_VNA_SHARED_H
#include <stddef.h>
#include <mpi.h>
/* FP32 ABI matches native Type_DS_VNA. Metadata is copied at initialization;
 * caller-owned atom/neighbor/image arrays may later be freed independently. */
int CpuVnaShared_Requested(void);
/* Call collectively and unconditionally before any flag-dependent Init. */
int CpuVnaShared_Agree(MPI_Comm comm);
int CpuVnaShared_Active(void);
void CpuVnaShared_Init(MPI_Comm comm,int atoms,const int*owners,const int*fan,
 int**natn,int**ncn,int**image4,const int*norb,int projectors,int scratch_fan,int scratch_norb);
float*CpuVnaShared_Row(int component,int atom,int neighbor,int orbital);
/* Epoch0: raw values and derivatives; epoch1: derivatives weighted by ene. */
void CpuVnaShared_Publish(int epoch);
/* All raw readers complete before owners modify derivative components. */
void CpuVnaShared_BeginWeighted(void);
void CpuVnaShared_BindAtom(float*****view,int scratch,int atom);
void CpuVnaShared_BindCenter(float*****view,int scratch,int center,int**reverse);
size_t CpuVnaShared_OwnedBytes(void);
void CpuVnaShared_Destroy(void);
#endif
