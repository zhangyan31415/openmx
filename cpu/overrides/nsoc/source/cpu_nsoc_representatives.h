#ifndef CPU_NSOC_REPRESENTATIVES_H
#define CPU_NSOC_REPRESENTATIVES_H
#include "cpu_nsoc_owner.h"
/* All ranks receive the same global representative descriptors, ordered by
 * native directed block ID. Ownership is the unchanged original atom owner. */
CpuNOBlock *CpuRSym_RepresentativeBlocks(int *count);
/* Per-rank values use CpuNOCompute's original-owner/global-ID order. */
void CpuRSym_RepresentativeExpand(const double *values,size_t entries,double ****field);
/* Validates correspondence between the solver's paired representation and
 * this projector's point-group generators; no coefficient communication. */
int CpuRSym_RepresentativeCertificate(int n,int maxn,const double *ko);
int CpuNSym_RepresentativeCertificate(int n,int maxn,const double *ko,
 int ng,int nat,const int *const *permutations,
 const int *types,const int *norb,const double *const *rotations);
int CpuNSym_RepresentativeTailCandidate(int n,int maxn,const double *ko,
 int ng,int nat,const int *const *permutations,const int *types,const int *norb,
 const double *const *rotations,int *cut_count,const int **cut_ids);
int CpuRSym_RepresentativeTailCandidate(int n,int maxn,const double *ko,
 int *cut_count,const int **cut_ids,double *amplification);
#endif
