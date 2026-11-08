#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* The planner's failure path is irrelevant to this standalone equivalence
 * test; production still uses its native collective abort. */
#define mpi_comm_level1 0
#define MPI_Abort(comm,code) ((void)(comm),(void)(code),abort())
#include "force3_tile_plan.h"
static int get(void *opaque,int h,int k){return ((int**)opaque)[h][k];}
int main(void)
{
  int a[]={0,1,6,7,7,8,13,14,21,22,56,62,63};
  int b[]={63,1,0,14,7,56,1,8,22,62,13,21,6};
  int z=0,*rows[]={a,b,&z},length[]={13,13,0};
  for(int bounded=0;bounded<2;bounded++) {
    size_t used0,used1,budget=bounded?3*sizeof(CpuF3TileIndex):1024*1024;
    int scan0,scan1;
    CpuF3TileIndex *old=CpuF3TilePrepare(3,length,rows,64,7,budget,&used0,&scan0);
    /* NULL proves that the indexed route never reads the raw pointer tree. */
    CpuF3TileIndex *now=CpuF3TilePrepareIndexed(3,length,NULL,64,7,budget,&used1,&scan1,get,rows);
    assert(used0==used1 && scan0==scan1);
    for(int h=0;h<3;h++) {
      assert(old[h].mode==now[h].mode);
      for(int t=0;t<10;t++) {
        int begin0,end0,begin1,end1,stop=(t+1)*7;if(stop>64)stop=64;
        CpuF3TileBounds(old+h,length[h],rows[h],t,stop,&begin0,&end0);
        CpuF3TileBoundsIndexed(now+h,length[h],NULL,t,stop,&begin1,&end1,get,rows,h);
        assert(begin0==begin1 && end0==end1);
        for(int k=begin0;k<end0;k++) {
          int q0=old[h].mode==1?old[h].order[k]:k;
          int q1=now[h].mode==1?now[h].order[k]:k;
          assert(q0==q1 && rows[h][q0]==get(rows,h,q1));
        }
      }
    }
    CpuF3TileFree(old,3);CpuF3TileFree(now,3);
  }
  puts("PASS Force3 raw/indexed plans preserve monotone, stable-bucket and bounded-rescan visit order");
  return 0;
}
