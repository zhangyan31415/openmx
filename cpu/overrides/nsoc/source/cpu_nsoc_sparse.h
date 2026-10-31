#ifndef OPENMX_CPU_NSOC_SPARSE_H
#define OPENMX_CPU_NSOC_SPARSE_H
#include <stdint.h>
#include <stddef.h>
#include "mpi.h"
typedef struct CpuNSparsePlan CpuNSparsePlan;
typedef struct {double right_pack,right_collective,right_multiply,left_pack,left_collective,left_multiply,total;} CpuNSparseTimings;
CpuNSparsePlan *CpuNSparse_Create(int n,int nb,int nr,int nc,int np_rows,int np_cols,int prow,int pcol,const int64_t*cp,const int*ri,const double*cv,MPI_Comm world,char*error,size_t error_size);
int CpuNSparse_Congruence(CpuNSparsePlan*p,double*A,int lda,double*work,int ldw,int profile,CpuNSparseTimings*timing,char*error,size_t error_size);
/* W=U*A, with the same N-by-N block-cyclic layout as the congruence.
 * Input/output must not alias. batch_rows=0 retains one component/collective;
 * 1..1024 groups whole components in a separate reusable bounded panel. */
int CpuNSparse_Left(CpuNSparsePlan*p,const double*A,int lda,double*work,int ldw,int batch_rows,int profile,CpuNSparseTimings*timing,char*error,size_t error_size);
size_t CpuNSparse_BackPanelBytes(const CpuNSparsePlan*p);
int CpuNSparse_BackBatches(const CpuNSparsePlan*p);
size_t CpuNSparse_PanelBytes(const CpuNSparsePlan*p);
size_t CpuNSparse_MetadataBytes(const CpuNSparsePlan*p);
int CpuNSparse_Components(const CpuNSparsePlan*p);
int CpuNSparse_MaxBlock(const CpuNSparsePlan*p);
void CpuNSparse_Destroy(CpuNSparsePlan*p);
#endif
