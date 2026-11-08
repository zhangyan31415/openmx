#ifndef OPENMX_FORCE3_TILE_PLAN_H
#define OPENMX_FORCE3_TILE_PLAN_H
#include <stdint.h>
#include <limits.h>
/* This metadata exists only for the currently owned centre atom. */
typedef struct {
  int mode; /* 0: checked monotone cursor, 1: stable tile buckets, 2: bounded rescan */
  int cursor;
  int *order;
  size_t *offset;
} CpuF3TileIndex;
static void CpuF3TileFail(const char *why) {
  fprintf(stderr,"<CPU_FORCE3_TILE> %s\n",why);fflush(stderr);
  MPI_Abort(mpi_comm_level1,92);abort();
}
static void *CpuF3Alloc(size_t bytes) {
  void *p=malloc(bytes?bytes:1);if(!p)CpuF3TileFail("allocation failed");return p;
}
static size_t CpuF3CheckedProduct(size_t a,size_t b) {
  if(b && a>SIZE_MAX/b)CpuF3TileFail("workspace size overflow");return a*b;
}
typedef int (*CpuF3IndexGetter)(void *context,int neighbor,int point);
static inline int CpuF3IndexRead(int **nc,CpuF3IndexGetter get,void *context,int h,int k){
  return get?get(context,h,k):nc[h][k];
}
static CpuF3TileIndex *CpuF3TilePrepareIndexed(int nh,const int *length,int **nc,int ng,int tile,
                                      size_t budget,size_t *used,int *rescans,
                                      CpuF3IndexGetter get,void *context) {
  size_t base=CpuF3CheckedProduct((size_t)nh,sizeof(CpuF3TileIndex));
  if(base>budget)CpuF3TileFail("centre-neighbour plan exceeds metadata budget");
  CpuF3TileIndex *p=CpuF3Alloc(base);memset(p,0,base);*used=base;*rescans=0;
  size_t nt=((size_t)ng+(size_t)tile-1)/(size_t)tile;
  int h;
  for(h=0;h<nh;h++) {
    int n=length[h],k,sorted=1;
    if(n<0)CpuF3TileFail("negative overlap length");
    for(k=0;k<n;k++) {
      if(CpuF3IndexRead(nc,get,context,h,k)<0 || CpuF3IndexRead(nc,get,context,h,k)>=ng)CpuF3TileFail("overlap index outside centre grid");
      if(k && CpuF3IndexRead(nc,get,context,h,k)<CpuF3IndexRead(nc,get,context,h,k-1))sorted=0;
    }
    if(sorted)continue;
    size_t bytes=CpuF3CheckedProduct((size_t)n,sizeof(int));
    size_t offsets=CpuF3CheckedProduct(nt+1,sizeof(size_t));
    /* count buckets in offset[1..nt], then fill in reverse using endpoint cursors.
       No second per-tile cursor allocation and no changes to native GList arrays. */
    if(bytes>SIZE_MAX-offsets || *used>budget || bytes+offsets>budget-*used) {
      p[h].mode=2;(*rescans)++;continue;
    }
    p[h].mode=1;p[h].order=CpuF3Alloc(bytes);p[h].offset=CpuF3Alloc(offsets);
    memset(p[h].offset,0,offsets);
    for(k=0;k<n;k++)p[h].offset[(size_t)CpuF3IndexRead(nc,get,context,h,k)/(size_t)tile+1]++;
    for(size_t t=1;t<=nt;t++)p[h].offset[t]+=p[h].offset[t-1];
    for(k=n;k-->0;) {
      size_t t=(size_t)CpuF3IndexRead(nc,get,context,h,k)/(size_t)tile;
      p[h].order[--p[h].offset[t+1]]=k;
    }
    /* In-place endpoint cursors now contain starts shifted by one bucket. */
    for(size_t t=0;t<nt;t++)p[h].offset[t]=p[h].offset[t+1];
    p[h].offset[nt]=(size_t)n;
    *used+=bytes+offsets;
  }
  return p;
}
static void CpuF3TileBoundsIndexed(CpuF3TileIndex *p,int n,const int *nc,int tile_no,int end,
                            int *begin,int *stop,CpuF3IndexGetter get,void *context,int h) {
  if(p->mode==0) {
    *begin=p->cursor;
    while(p->cursor<n && (get?get(context,h,p->cursor):nc[p->cursor])<end)p->cursor++;
    *stop=p->cursor;
  } else if(p->mode==1) {
    *begin=(int)p->offset[tile_no];*stop=(int)p->offset[tile_no+1];
  } else {*begin=0;*stop=n;}
}
/* Native array API retained for the SOC/legacy paths and existing tests. */
static CpuF3TileIndex *CpuF3TilePrepare(int nh,const int *length,int **nc,int ng,int tile,
                                      size_t budget,size_t *used,int *rescans){
  return CpuF3TilePrepareIndexed(nh,length,nc,ng,tile,budget,used,rescans,NULL,NULL);
}
static void CpuF3TileBounds(CpuF3TileIndex *p,int n,const int *nc,int tile_no,int end,
                            int *begin,int *stop){
  CpuF3TileBoundsIndexed(p,n,nc,tile_no,end,begin,stop,NULL,NULL,0);
}
static void CpuF3TileFree(CpuF3TileIndex *p,int nh) {
  for(int h=0;h<nh;h++){free(p[h].order);free(p[h].offset);}free(p);
}
#endif
