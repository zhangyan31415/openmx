#ifndef OPENMX_CPU_OVERLAP_INDEX_H
#define OPENMX_CPU_OVERLAP_INDEX_H
#include <stdio.h>
#include <stdlib.h>
#include "cpu_overlap_runs.h"
#include "cpu_grid_shared.h"
/* Include openmx_common.h first. Raw native pointer trees remain the fallback
 * for unsupported modes and for pairs whose runs would use more memory. */
extern CpuOverlapPair **CpuOverlapPairs;
extern uint64_t CpuOverlapGeneration;
int CpuOverlapEligible(void);
int CpuSOCMemoryEligible(void);
void CpuOverlapInit(void);
void CpuOverlapStore(int mc,int h,int count,const int *nc,const int *nh);
void CpuOverlapReport(void);
void CpuOverlapDestroy(void);
static __thread CpuOverlapCursor cpu_overlap_cursor;
static inline int CpuOverlapIndex(int axis,int mc,int h,int nog)
{
  if(CpuOverlapPairs) {
    const CpuOverlapPair *p=&CpuOverlapPairs[mc][h];
    if(p->runs) {
      const CpuOverlapRun *r=CpuOverlapRunFind(p,CpuOverlapGeneration,nog,
                                             &cpu_overlap_cursor);
      if(!r) { fprintf(stderr,"CPU_OVERLAP_INDEX invalid access %d %d %d\n",mc,h,nog);abort(); }
      return (axis==0?r->nc:r->nh)+(nog-r->begin);
    }
  }
  return axis==0?GListTAtoms1[mc][h][nog]:GListTAtoms2[mc][h][nog];
}
static inline int CpuOverlapIndex1(int mc,int h,int nog){return CpuOverlapIndex(0,mc,h,nog);}
static inline int CpuOverlapIndex2(int mc,int h,int nog){return CpuOverlapIndex(1,mc,h,nog);}
/* Same-node orbital values are already shared. Their row address can be
 * derived from the existing owner-grid index instead of saving another pointer
 * for every overlap occurrence. Off-node values retain native private rows. */
static inline Type_Orbs_Grid *CpuOverlapRow(int mc,int h,int nog)
{
  int atom=natn[M2G[mc]][h];
  if(CpuOverlapPairs && CpuGridShared_IsAtom(atom))
    return CpuGridShared_ReadRow(atom,CpuOverlapIndex2(mc,h,nog));
  return Orbs_Grid_FNAN[mc][h][nog];
}
void CpuGridLifetime_Release(void);
#endif
