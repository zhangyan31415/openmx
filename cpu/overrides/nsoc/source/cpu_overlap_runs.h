#ifndef OPENMX_CPU_OVERLAP_RUNS_H
#define OPENMX_CPU_OVERLAP_RUNS_H
#include <stddef.h>
#include <stdint.h>
#include <limits.h>

/* An exact encoding of the existing Nog ordering. A run never reorders a
 * quadrature point: both original local indices advance by exactly one. */
typedef struct { int begin, nc, nh, length; } CpuOverlapRun;
typedef struct {
  int count, nruns;
  CpuOverlapRun *runs;
} CpuOverlapPair;
typedef struct {
  const CpuOverlapPair *pair;
  uint64_t generation;
  int run;
} CpuOverlapCursor;

static inline int CpuOverlapRunCount(int n,const int *nc,const int *nh)
{
  int runs=n>0?1:0;
  for(int k=1;k<n;k++)
    if(nc[k-1]==INT_MAX || nh[k-1]==INT_MAX ||
       nc[k]!=nc[k-1]+1 || nh[k]!=nh[k-1]+1) ++runs;
  return runs;
}
static inline void CpuOverlapRunEncode(int n,const int *nc,const int *nh,
                                       CpuOverlapRun *out)
{
  int r=-1;
  for(int k=0;k<n;k++) {
    if(k==0 || nc[k-1]==INT_MAX || nh[k-1]==INT_MAX ||
       nc[k]!=nc[k-1]+1 || nh[k]!=nh[k-1]+1) {
      ++r;out[r]=(CpuOverlapRun){k,nc[k],nh[k],1};
    } else ++out[r].length;
  }
}
/* The caller owns the cursor; production uses one cursor per OS thread.
 * Random access and revisiting a previous pair are exact, without a linear
 * scan through all preceding runs. Generation protects reinitialization. */
static inline const CpuOverlapRun *CpuOverlapRunFind(const CpuOverlapPair *p,
    uint64_t generation,int nog,CpuOverlapCursor *cursor)
{
  int r;
  if(!p || !p->runs || nog<0 || nog>=p->count || p->nruns<1) return NULL;
  if(cursor->pair==p && cursor->generation==generation) {
    r=cursor->run;
    if(r>=0 && r<p->nruns) {
      const CpuOverlapRun *q=p->runs+r;
      if(nog>=q->begin && nog-q->begin<q->length) return q;
      if(r+1<p->nruns && nog==q->begin+q->length) {
        cursor->run=r+1;return p->runs+r+1;
      }
    }
  }
  int lo=0,hi=p->nruns;
  while(lo+1<hi) {
    int mid=lo+(hi-lo)/2;
    if(p->runs[mid].begin<=nog)lo=mid;else hi=mid;
  }
  cursor->pair=p;cursor->generation=generation;cursor->run=lo;
  return p->runs+lo;
}
#endif
