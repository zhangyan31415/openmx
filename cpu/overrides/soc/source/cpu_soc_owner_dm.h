#ifndef OPENMX_CPU_SOC_OWNER_DM_H
#define OPENMX_CPU_SOC_OWNER_DM_H
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "owner_sparse.h"
#include "cpu_dm_compact_trim.h"
#include "cpu_dm_balance_config.h"
/* Caller provides mpi.h and cpu_dm_blas3.v2.h. Metadata is identical on ranks.
 * Blocks retain native global output offsets; only owned values are retained.
 * Temporary component output + BLAS packing/scratch is bounded by tile_bytes.
 */
typedef struct {
  int row, rows, col, cols, offset, owner;
} CpuSocOwnerBlock;
typedef struct {
  double *values;
  int *global_offset, *length;
  int blocks, entries, components;
  size_t cache_bytes, peak_tile_bytes;
} CpuSocOwnerCache;
static void cpu_soc_owner_free(CpuSocOwnerCache *c)
{
  free(c->values); free(c->global_offset); free(c->length);
  memset(c,0,sizeof(*c));
}
/* The kernel allocates G[maxm*maxc] with independent cross-group maxima.
 * Use the validated global maximum row count so tall and wide groups cannot
 * defeat the sum-of-blocks bound. It also bounds coalesced BLAS packing.
 * Block metadata for the tile is included; existing eigenvectors are excluded. */
static size_t cpu_soc_owner_block_bytes(const CpuSocOwnerBlock *b,int comps,int band_tile,int max_rows)
{
  return (size_t)comps*b->rows*b->cols*sizeof(double)+
    (size_t)16*((size_t)band_tile*((size_t)b->rows+b->cols)+(size_t)max_rows*b->cols)+5*sizeof(int);
}
static int cpu_soc_owner_collective_error(int rc,MPI_Comm comm)
{
  int bad=rc?1:0, any=0;
  MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,comm);
  return any ? (rc?rc:-8) : 0;
}
static int cpu_soc_owner_build(CpuSocOwnerCache *cache,int n,int bands,
 const void *vectors,const double *weights,int count,const CpuSocOwnerBlock *blocks,
 int global_entries,int edm,size_t tile_bytes,MPI_Comm comm)
{
  int rank,np,b,owner,rc=0,*order=NULL,*start=NULL,*next=NULL,*tile=NULL;
  int owned_blocks=0,owned_entries=0,max_tile_blocks=0,max_tile_entries=0;
  int comps=edm?4:6,band_tile=256,max_rows=0;
  double *partial=NULL;
  const int profile=CpuDmBalanceFlag("OPENMX_CPU_DM_BALANCE_PROFILE")==1;
  CpuDMCompactProfile total={0};
  double comm_seconds=0.,comm_start=0.;
  size_t owned_bytes=0,max_kernel_extra=0,tile_budget=tile_bytes/2;
  MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&np);
  cpu_soc_owner_free(cache);
  if(n<0||bands<0||count<0||global_entries<0||(count&&!blocks)||
     (bands&&(!vectors||!weights))||!(edm==0||edm==1)||tile_bytes<2)rc=-1;
  if(cpu_soc_owner_collective_error(rc,comm))return -1;
  start=(int*)calloc((size_t)np+1,sizeof(int));
  next=(int*)calloc((size_t)np,sizeof(int));
  order=(int*)malloc((size_t)(count?count:1)*sizeof(int));
  if(!start||!next||!order)rc=-2;
  if(cpu_soc_owner_collective_error(rc,comm)){rc=-2;goto done;}
  for(b=0;b<count;b++){
    const CpuSocOwnerBlock *x=blocks+b;
    int64_t len=(int64_t)x->rows*x->cols;
    if(x->rows<=0||x->cols<=0||x->row<0||x->col<0||x->offset<0||
       (int64_t)x->row+x->rows>n||(int64_t)x->col+x->cols>n||
       (int64_t)x->offset+len>global_entries||x->owner<0||x->owner>=np||
       len>INT_MAX/comps){rc=-3;break;}
    if(x->rows>max_rows)max_rows=x->rows;
    start[x->owner+1]++;
    if(x->owner==rank){
      owned_blocks++;
      if((int64_t)owned_entries+len>INT_MAX){rc=-3;break;}
      owned_entries+=(int)len;
    }
  }
  if(cpu_soc_owner_collective_error(rc,comm)){rc=-3;goto done;}
  for(b=0;b<count;b++)if(cpu_soc_owner_block_bytes(blocks+b,comps,band_tile,max_rows)>tile_budget){rc=-3;break;}
  if(cpu_soc_owner_collective_error(rc,comm)){rc=-3;goto done;}
  for(owner=0;owner<np;owner++){start[owner+1]+=start[owner];next[owner]=start[owner];}
  for(b=0;b<count;b++)order[next[blocks[b].owner]++]=b;
  /* Plan identical owner-ordered tiles on every rank, including empty owners. */
  for(owner=0;owner<np;owner++){
    int pos=start[owner];
    while(pos<start[owner+1]){
      int nt=0,ne=0;size_t cost=0;
      while(pos<start[owner+1]){
        const CpuSocOwnerBlock *x=blocks+order[pos];
        size_t add=cpu_soc_owner_block_bytes(x,comps,band_tile,max_rows);
        int len=x->rows*x->cols;
        if(nt&&(add>tile_budget-cost||len>INT_MAX/comps-ne))break;
        cost+=add;ne+=len;nt++;pos++;
      }
      if(nt>max_tile_blocks)max_tile_blocks=nt;
      if(ne>max_tile_entries)max_tile_entries=ne;
      if(cost-(size_t)comps*ne*sizeof(double)-(size_t)5*nt*sizeof(int)>max_kernel_extra)
        max_kernel_extra=cost-(size_t)comps*ne*sizeof(double)-(size_t)5*nt*sizeof(int);
    }
  }
  cache->peak_tile_bytes=(size_t)comps*max_tile_entries*sizeof(double)+max_kernel_extra+(size_t)5*max_tile_blocks*sizeof(int);
  /* The output and metadata buffers are reusable maxima across all tiles,
   * so bound those independently from the kernel's maximum scratch demand. */
  if(cache->peak_tile_bytes>tile_bytes)rc=-3;
  if(cpu_soc_owner_collective_error(rc,comm)){rc=-3;goto done;}
  owned_bytes=(size_t)comps*owned_entries*sizeof(double);
  cache->values=(double*)malloc(owned_bytes?owned_bytes:sizeof(double));
  cache->global_offset=(int*)malloc((size_t)(owned_blocks?owned_blocks:1)*sizeof(int));
  cache->length=(int*)malloc((size_t)(owned_blocks?owned_blocks:1)*sizeof(int));
  partial=(double*)malloc((size_t)(max_tile_entries?max_tile_entries:1)*comps*sizeof(double));
  tile=(int*)malloc((size_t)(max_tile_blocks?max_tile_blocks:1)*5*sizeof(int));
  if(!cache->values||!cache->global_offset||!cache->length||!partial||!tile)rc=-2;
  if(cpu_soc_owner_collective_error(rc,comm)){rc=-2;goto done;}
  cache->blocks=owned_blocks;cache->entries=owned_entries;cache->components=comps;
  cache->cache_bytes=owned_bytes+(size_t)2*owned_blocks*sizeof(int);
  {
    int at=0;
    for(b=start[rank];b<start[rank+1];b++){
      const CpuSocOwnerBlock *x=blocks+order[b];
      cache->global_offset[at]=x->offset;cache->length[at]=x->rows*x->cols;at++;
    }
  }
  for(owner=0;owner<np;owner++){
    int pos=start[owner],owner_at=0;
    while(pos<start[owner+1]){
      int nt=0,ne=0,c;size_t cost=0;
      while(pos<start[owner+1]){
        const CpuSocOwnerBlock *x=blocks+order[pos];
        size_t add=cpu_soc_owner_block_bytes(x,comps,band_tile,max_rows);
        int len=x->rows*x->cols;
        if(nt&&(add>tile_budget-cost||len>INT_MAX/comps-ne))break;
        tile[5*nt]=x->row;tile[5*nt+1]=x->rows;
        tile[5*nt+2]=x->col;tile[5*nt+3]=x->cols;tile[5*nt+4]=ne;
        cost+=add;ne+=len;nt++;pos++;
      }
      rc=CpuDMCoefficient_Enabled()?cpu_dm_blocks_compact(n,bands,vectors,weights,nt,tile,ne,edm,band_tile,partial):cpu_dm_blocks(n,bands,vectors,weights,nt,tile,ne,edm,band_tile,partial);
      if(profile){
        CpuDMCompactProfile p;CpuDMCompactGetProfile(&p);
        total.pack_calls+=p.pack_calls;total.gemm_calls+=p.gemm_calls;
        total.kernel_seconds+=p.kernel_seconds;total.pack_seconds+=p.pack_seconds;
        total.gemm_seconds+=p.gemm_seconds;total.scatter_seconds+=p.scatter_seconds;
        if(p.scratch_bytes>total.scratch_bytes)total.scratch_bytes=p.scratch_bytes;
        comm_start=MPI_Wtime();
      }
      if(cpu_soc_owner_collective_error(rc,comm)){rc=-4;goto done;}
      if(rank==owner)MPI_Reduce(MPI_IN_PLACE,partial,comps*ne,MPI_DOUBLE,MPI_SUM,owner,comm);
      else MPI_Reduce(partial,partial,comps*ne,MPI_DOUBLE,MPI_SUM,owner,comm);
      if(profile)comm_seconds+=MPI_Wtime()-comm_start;
      if(rank==owner){
        for(c=0;c<comps;c++)memcpy(cache->values+(size_t)c*owned_entries+owner_at,
                                   partial+(size_t)c*ne,(size_t)ne*sizeof(double));
        owner_at+=ne;
      }
    }
  }
 done:
  if(!rc&&profile)printf("DM_BALANCE_PROFILE rank=%d edm=%d pack_calls=%zu gemm_calls=%zu scratch_bytes=%zu kernel_seconds=%.9g pack_seconds=%.9g gemm_seconds=%.9g scatter_seconds=%.9g collective_seconds=%.9g\n",
    rank,edm,total.pack_calls,total.gemm_calls,total.scratch_bytes,total.kernel_seconds,
    total.pack_seconds,total.gemm_seconds,total.scatter_seconds,comm_seconds);
  free(start);free(next);free(order);free(tile);free(partial);
  if(rc)cpu_soc_owner_free(cache);
  return rc;
}
static int cpu_soc_owner_scatter(const CpuSocOwnerCache *c,int component,int global_entries,double *out)
{
  int b,at=0;
  if(component<0||component>=c->components||global_entries<0||!out||!c->values)return -1;
  if(!CpuOwnerSparseRequested())memset(out,0,(size_t)global_entries*sizeof(double));
  for(b=0;b<c->blocks;b++){
    if(c->global_offset[b]<0||(int64_t)c->global_offset[b]+c->length[b]>global_entries)return -1;
    memcpy(out+c->global_offset[b],c->values+(size_t)component*c->entries+at,(size_t)c->length[b]*sizeof(double));
    at+=c->length[b];
  }
  return at==c->entries?0:-1;
}
static int cpu_soc_owner_compare(const CpuSocOwnerCache *c,const double *reference,const double *candidate,double *error,double *scale)
{
  int b,j;*error=0;*scale=0;
  for(b=0;b<c->blocks;b++)for(j=0;j<c->length[b];j++){
    int at=c->global_offset[b]+j;double x=reference[at],y=candidate[at],e=fabs(x-y);
    if(!isfinite(x)||!isfinite(y))return -1;
    if(e>*error)*error=e;if(fabs(x)>*scale)*scale=fabs(x);
  }
  return 0;
}
#endif
