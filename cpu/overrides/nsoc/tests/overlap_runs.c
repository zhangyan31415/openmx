#include <assert.h>
#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "cpu_overlap_runs.h"

enum { N=150013, THREADS=8 };
static int nc[N],nh[N];
static CpuOverlapPair pair;
static __thread CpuOverlapCursor cursor;

static void check(int k,uint64_t generation)
{
  const CpuOverlapRun *r=CpuOverlapRunFind(&pair,generation,k,&cursor);
  assert(r && k>=r->begin && k-r->begin<r->length);
  assert(r->nc+k-r->begin==nc[k]);
  assert(r->nh+k-r->begin==nh[k]);
}
static void *worker(void *opaque)
{
  uintptr_t id=(uintptr_t)opaque;
  /* Interleave sequential, reverse, same-point and random requests. The
   * cursor must stay independent when different threads read different runs. */
  for(int k=0;k<N;k++)check(k,1);
  for(int k=N;k-->0;)check(k,1);
  uint32_t state=(uint32_t)id+17;
  for(int k=0;k<200000;k++) {
    state=1664525u*state+1013904223u;
    int q=(int)(state%N);check(q,1);check(q,1);
  }
  return NULL;
}
int main(void)
{
  int a=0,b=71381;
  assert(CpuOverlapRunCount(0,NULL,NULL)==0);
  assert(CpuOverlapRunCount(1,&a,&b)==1);
  for(int k=0;k<N;k++) {
    /* Runs cross the old uint16 boundary; irregular jumps, reverse ranges,
     * duplicate values and single-point runs must all survive exactly. */
    int block=k/127;
    nc[k]=71000+(block%3==0?-1:1)*(k%127)+block*137;
    nh[k]=32000+(k%127)+block*139;
    if(k%193==0)nc[k]=3;
  }
  pair.count=N;pair.nruns=CpuOverlapRunCount(N,nc,nh);
  pair.runs=malloc((size_t)pair.nruns*sizeof(*pair.runs));assert(pair.runs);
  CpuOverlapRunEncode(N,nc,nh,pair.runs);
  assert(CpuOverlapRunFind(&pair,1,-1,&cursor)==NULL);
  assert(CpuOverlapRunFind(&pair,1,N,&cursor)==NULL);
  pthread_t threads[THREADS];
  for(uintptr_t t=0;t<THREADS;t++)assert(pthread_create(threads+t,NULL,worker,(void*)t)==0);
  for(int t=0;t<THREADS;t++)assert(pthread_join(threads[t],NULL)==0);
  /* Same address, a different encoding: generation must invalidate cursors. */
  check(1000,1);
  for(int k=0;k<N;k++){nc[k]=k;nh[k]=k+100;}
  pair.nruns=CpuOverlapRunCount(N,nc,nh);assert(pair.nruns==1);
  CpuOverlapRunEncode(N,nc,nh,pair.runs);
  check(1000,2);check(N-1,2);
  int edge[3]={INT_MAX-1,INT_MAX,0},other[3]={4,5,6};
  CpuOverlapRun out[3];assert(CpuOverlapRunCount(3,edge,other)==2);
  CpuOverlapRunEncode(3,edge,other,out);assert(out[0].length==2 && out[1].begin==2);
  printf("PASS exact overlap run boundaries, 16-bit overflow, random/reverse access, %d concurrent cursors and reinitialization\n",THREADS);
  free(pair.runs);return 0;
}
