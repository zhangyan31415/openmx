#ifndef CPU_SPARSE_LIFETIME_H
#define CPU_SPARSE_LIFETIME_H
#include <stddef.h>
#include <stdlib.h>

/* Retire whole native history slots after their last SCF reader. The outer
 * array and the first `keep` slots remain owned by the native caller. Null
 * slots make later cleanup idempotent; this does not shorten SCF history. */
static size_t CpuSparseHistoryRelease(double ******tree, int slots, int keep,
    int channels, int atoms, const int *map, const int *species,
    const int *norb, const int *fan)
{
  size_t rows=0;
  if(!tree)return 0;
  for(int s=keep;s<slots;s++){
    if(!tree[s])continue;
    for(int k=0;k<channels;k++){
      for(int a=0;a<=atoms;a++){
        int g=a?map[a]:0, n=a?norb[species[g]]:1, f=a?fan[g]:0;
        for(int h=0;h<=f;h++){
          for(int i=0;i<n;i++){free(tree[s][k][a][h][i]);rows++;}
          free(tree[s][k][a][h]);
        }
        free(tree[s][k][a]);
      }
      free(tree[s][k]);
    }
    free(tree[s]);tree[s]=NULL;
  }
  return rows;
}
#endif
