#include "cpu_nsoc_owner.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <mkl_cblas.h>
#ifndef CPU_NO_STAR_CACHE_BYTES
#define CPU_NO_STAR_CACHE_BYTES (64u*1024u*1024u)
#endif
struct CpuNOPlan{MPI_Comm comm;int rank,np,nao,nb,maxcols;CpuNOBlock*b;int*order,*start;size_t cap,max_tile,owned,weighted_cap,provider_peak;double *partial,*weighted,*d,*ed;};
static int any(MPI_Comm c,int bad){int x;MPI_Allreduce(&bad,&x,1,MPI_INT,MPI_MAX,c);return x;}
static uint64_t hashword(uint64_t h,int x){return(h^(uint32_t)x)*UINT64_C(1099511628211);}
void CpuNODestroy(CpuNOPlan*p){if(!p)return;free(p->b);free(p->order);free(p->start);free(p->partial);free(p->weighted);free(p->d);free(p->ed);free(p);}
CpuNOPlan*CpuNOCreate(MPI_Comm comm,int nao,int nb,const CpuNOBlock*b,size_t tile_bytes){
 int rank,np,bad=0,m[2]={nao,nb},lo[2],hi[2];unsigned long long cap=tile_bytes,cl,ch;
 MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&np);MPI_Allreduce(m,lo,2,MPI_INT,MPI_MIN,comm);MPI_Allreduce(m,hi,2,MPI_INT,MPI_MAX,comm);MPI_Allreduce(&cap,&cl,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,comm);MPI_Allreduce(&cap,&ch,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm);
 if(nao<=0||nb<0||(nb&&!b)||lo[0]!=hi[0]||lo[1]!=hi[1]||cl!=ch||tile_bytes<2*sizeof(double)||tile_bytes/sizeof(double)>INT_MAX)bad=1;
 if(any(comm,bad))return NULL;
 CpuNOPlan*p=calloc(1,sizeof(*p));if(any(comm,!p)){free(p);return NULL;}p->comm=comm;p->rank=rank;p->np=np;p->nao=nao;p->nb=nb;p->cap=tile_bytes/(2*sizeof(double));
 p->b=malloc((size_t)(nb?nb:1)*sizeof(*b));p->order=malloc((size_t)(nb?nb:1)*sizeof(int));p->start=calloc((size_t)np+1,sizeof(int));
 if(any(comm,!p->b||!p->order||!p->start)){CpuNODestroy(p);return NULL;}
 uint64_t hash=UINT64_C(1469598103934665603),hl,hh;
 for(int k=0;k<nb;k++){
  const CpuNOBlock*x=b+k;int64_t nr=(int64_t)x->row+x->rows,nc=(int64_t)x->col+x->cols;
  if(x->row<0||x->col<0||x->rows<=0||x->cols<=0||nr>nao||nc>nao||x->owner<0||x->owner>=np){bad=1;continue;}
  size_t len=(size_t)x->rows*x->cols;if(len>p->cap){bad=1;continue;}
  hash=hashword(hash,x->row);hash=hashword(hash,x->rows);hash=hashword(hash,x->col);hash=hashword(hash,x->cols);hash=hashword(hash,x->owner);
  p->b[k]=*x;p->start[x->owner+1]++;if(x->cols>p->maxcols)p->maxcols=x->cols;
  if(x->owner==rank){if(p->owned>SIZE_MAX/(2*sizeof(double))-len)bad=1;else p->owned+=len;}
 }
 if(any(comm,bad)){CpuNODestroy(p);return NULL;}
 MPI_Allreduce(&hash,&hl,1,MPI_UINT64_T,MPI_MIN,comm);MPI_Allreduce(&hash,&hh,1,MPI_UINT64_T,MPI_MAX,comm);if(hl!=hh){CpuNODestroy(p);return NULL;}
 for(int k=0;k<np;k++)p->start[k+1]+=p->start[k];int*next=malloc((size_t)np*sizeof(int));if(any(comm,!next)){free(next);CpuNODestroy(p);return NULL;}memcpy(next,p->start,(size_t)np*sizeof(int));
 for(int k=0;k<nb;k++)p->order[next[b[k].owner]++]=k;free(next);
 for(int owner=0;owner<np;owner++){int pos=p->start[owner];while(pos<p->start[owner+1]){size_t ne=0;while(pos<p->start[owner+1]){CpuNOBlock*x=p->b+p->order[pos];size_t len=(size_t)x->rows*x->cols;if(ne&&len>p->cap-ne)break;ne+=len;pos++;}if(ne>p->max_tile)p->max_tile=ne;}}
 p->partial=malloc((p->max_tile?2*p->max_tile:1)*sizeof(double));p->d=malloc((p->owned?p->owned:1)*sizeof(double));p->ed=malloc((p->owned?p->owned:1)*sizeof(double));
 if(any(comm,!p->partial||!p->d||!p->ed)){CpuNODestroy(p);return NULL;}return p;
}
int CpuNOCompute(CpuNOPlan*p,int nk,const double*c,const double*f,const double*energy){
 if(!p)return -1;int bad=nk<0||(nk&&(!c||!f||!energy));
 if(!bad)for(int k=0;k<nk;k++)if(!isfinite(f[k])||!isfinite(energy[k]))bad=1;
 if(!bad&&nk&&(size_t)p->nao>SIZE_MAX/sizeof(double)/(size_t)nk)bad=1;
 if(any(p->comm,bad))return -1;
 size_t need=(size_t)p->maxcols*(size_t)nk;
 if(need>p->weighted_cap){double*w=realloc(p->weighted,need*sizeof(double));if(!w)bad=1;else{p->weighted=w;p->weighted_cap=need;}}
 if(any(p->comm,bad))return -2;
 for(int owner=0;owner<p->np;owner++){
  int pos=p->start[owner];size_t owner_at=0;
  while(pos<p->start[owner+1]){
   int first=pos;size_t ne=0;
   while(pos<p->start[owner+1]){CpuNOBlock*x=p->b+p->order[pos];size_t len=(size_t)x->rows*x->cols;if(ne&&len>p->cap-ne)break;ne+=len;pos++;}
   size_t at=0;
   if(!nk)memset(p->partial,0,2*ne*sizeof(double));
   else for(int bi=first;bi<pos;bi++){
    CpuNOBlock*x=p->b+p->order[bi];const double*a=c+(size_t)x->row*nk,*b=c+(size_t)x->col*nk;
    for(int j=0;j<x->cols;j++)for(int k=0;k<nk;k++)p->weighted[(size_t)j*nk+k]=b[(size_t)j*nk+k]*f[k];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,x->rows,x->cols,nk,1.,a,nk,p->weighted,nk,0.,p->partial+at,x->cols);
    for(int j=0;j<x->cols;j++)for(int k=0;k<nk;k++)p->weighted[(size_t)j*nk+k]=b[(size_t)j*nk+k]*f[k]*energy[k];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,x->rows,x->cols,nk,1.,a,nk,p->weighted,nk,0.,p->partial+ne+at,x->cols);
    at+=(size_t)x->rows*x->cols;
   }
   int count=(int)(2*ne),rc;
   if(p->rank==owner)rc=MPI_Reduce(MPI_IN_PLACE,p->partial,count,MPI_DOUBLE,MPI_SUM,owner,p->comm);
   else rc=MPI_Reduce(p->partial,NULL,count,MPI_DOUBLE,MPI_SUM,owner,p->comm);
   if(rc!=MPI_SUCCESS)MPI_Abort(p->comm,rc);
   if(p->rank==owner){memcpy(p->d+owner_at,p->partial,ne*sizeof(double));memcpy(p->ed+owner_at,p->partial+ne,ne*sizeof(double));owner_at+=ne;}
  }
 }
 return 0;
}
typedef struct {int row,rows;double*data;} CpuNOCachedRows;
typedef struct {CpuNOCachedRows*entry;int count,capacity,center;size_t bytes,peak;int nk;const CpuNORowProvider*provider;MPI_Comm comm;} CpuNOStar;
static void CpuNOStarClear(CpuNOStar*s){for(int i=0;i<s->count;i++)free(s->entry[i].data);s->count=0;s->bytes=0;}
static const double*CpuNOStarRows(CpuNOStar*s,int row,int rows){
 for(int i=0;i<s->count;i++)if(s->entry[i].row==row&&s->entry[i].rows==rows)return s->entry[i].data;
 size_t bytes=(size_t)rows*s->nk*sizeof(double);
 /* Keep the center stable while fetching a neighbor. Only inactive neighbor
  * cache entries are evicted; the kernel never keeps a pointer past a block. */
 if(s->bytes+bytes>CPU_NO_STAR_CACHE_BYTES){
  int at=0;for(int i=0;i<s->count;i++)if(s->entry[i].row==s->center)s->entry[at++]=s->entry[i];else{free(s->entry[i].data);s->bytes-=(size_t)s->entry[i].rows*s->nk*sizeof(double);}s->count=at;
 }
 if(s->count==s->capacity){int capacity=s->capacity?s->capacity*2:16;CpuNOCachedRows*p=realloc(s->entry,(size_t)capacity*sizeof(*p));if(!p){MPI_Abort(s->comm,98);abort();}s->entry=p;s->capacity=capacity;}
 double*data=malloc(bytes?bytes:sizeof(double));if(!data){MPI_Abort(s->comm,98);abort();}
 if(s->provider->rows(s->provider->context,row,rows,data,s->nk)){MPI_Abort(s->comm,98);abort();}
 s->entry[s->count++]=(CpuNOCachedRows){row,rows,data};s->bytes+=bytes;if(s->bytes>s->peak)s->peak=s->bytes;return data;
}
int CpuNOComputeProvided(CpuNOPlan*p,int nk,const CpuNORowProvider*provider,const double*f,const double*energy){
 if(!p)return -1;int bad=nk<0||!provider||!provider->rows||(nk&&(!f||!energy));
 if(!bad)for(int k=0;k<nk;k++)if(!isfinite(f[k])||!isfinite(energy[k]))bad=1;
 if(!bad&&nk&&(size_t)p->nao>SIZE_MAX/sizeof(double)/(size_t)nk)bad=1;
 if(any(p->comm,bad))return -1;
 size_t need=(size_t)p->maxcols*nk;if(need>p->weighted_cap){double*w=realloc(p->weighted,need*sizeof(double));if(!w)bad=1;else{p->weighted=w;p->weighted_cap=need;}}
 if(any(p->comm,bad))return -2;
 CpuNOStar star={0};star.center=-1;star.nk=nk;star.provider=provider;star.comm=p->comm;
 for(int owner=0;owner<p->np;owner++){
  int pos=p->start[owner];size_t owner_at=0;
  while(pos<p->start[owner+1]){
   int first=pos;size_t ne=0;while(pos<p->start[owner+1]){CpuNOBlock*x=p->b+p->order[pos];size_t len=(size_t)x->rows*x->cols;if(ne&&len>p->cap-ne)break;ne+=len;pos++;}
   size_t at=0;if(!nk)memset(p->partial,0,2*ne*sizeof(double));
   else for(int bi=first;bi<pos;bi++){
    CpuNOBlock*x=p->b+p->order[bi];if(star.center!=x->row){CpuNOStarClear(&star);star.center=x->row;}
    const double*a=CpuNOStarRows(&star,x->row,x->rows),*b=CpuNOStarRows(&star,x->col,x->cols);
    for(int j=0;j<x->cols;j++)for(int k=0;k<nk;k++)p->weighted[(size_t)j*nk+k]=b[(size_t)j*nk+k]*f[k];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,x->rows,x->cols,nk,1.,a,nk,p->weighted,nk,0.,p->partial+at,x->cols);
    for(int j=0;j<x->cols;j++)for(int k=0;k<nk;k++)p->weighted[(size_t)j*nk+k]=b[(size_t)j*nk+k]*f[k]*energy[k];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,x->rows,x->cols,nk,1.,a,nk,p->weighted,nk,0.,p->partial+ne+at,x->cols);
    at+=(size_t)x->rows*x->cols;
   }
   int count=(int)(2*ne),rc;
   if(p->rank==owner)rc=MPI_Reduce(MPI_IN_PLACE,p->partial,count,MPI_DOUBLE,MPI_SUM,owner,p->comm);
   else rc=MPI_Reduce(p->partial,NULL,count,MPI_DOUBLE,MPI_SUM,owner,p->comm);
   if(rc!=MPI_SUCCESS)MPI_Abort(p->comm,rc);
   if(p->rank==owner){memcpy(p->d+owner_at,p->partial,ne*sizeof(double));memcpy(p->ed+owner_at,p->partial+ne,ne*sizeof(double));owner_at+=ne;}
  }
 }
 p->provider_peak=star.peak;CpuNOStarClear(&star);free(star.entry);return 0;
}
size_t CpuNOEntries(const CpuNOPlan*p){return p?p->owned:0;}
size_t CpuNOTemporaryBytes(const CpuNOPlan*p){return p?(2*p->max_tile+p->weighted_cap)*sizeof(double)+p->provider_peak:0;}
const double*CpuNODensity(const CpuNOPlan*p){return p?p->d:NULL;}
const double*CpuNOEnergyDensity(const CpuNOPlan*p){return p?p->ed:NULL;}
