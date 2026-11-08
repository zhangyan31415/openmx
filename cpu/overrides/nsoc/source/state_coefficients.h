#ifndef OPENMX_CPU_STATE_COEFFICIENTS_H
#define OPENMX_CPU_STATE_COEFFICIENTS_H
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sector_native_pack.h"
typedef struct {
 int n,nk,first,last,ns,ends[48],starts[48],count[48],position_start[49];
 size_t offset[48],used;int*rowptr,*columns,*positions,*state_sector,*state_slot,*source_slot;
 double*values,*data,*reference;double audit_error;
 const int64_t*source_cp;
} CpuStateCoefficients;
/* Retain sparse U in row order, with duplicate CSC terms collapsed in their
 * original order. Rows then use exactly the previous sparse-product sum. */
static int CpuStateInit(CpuStateCoefficients*c,int n,const int64_t*cp,const int*ri,const double*cv){
 if(c->n)return c->n==n&&c->source_cp==cp?0:-1;
 if(n<1||!cp||!ri||!cv||cp[n]<0||cp[n]>INT_MAX)return -1;
 int zmax=(int)cp[n];int*count=calloc((size_t)n+1,sizeof(int)),*cursor=calloc((size_t)n,sizeof(int));
 c->rowptr=calloc((size_t)n+1,sizeof(int));c->columns=malloc((size_t)(zmax?zmax:1)*sizeof(int));c->values=malloc((size_t)(zmax?zmax:1)*sizeof(double));
 if(!count||!cursor||!c->rowptr||!c->columns||!c->values){free(count);free(cursor);return -1;}
 for(int64_t z=0;z<cp[n];z++){if(ri[z]<0||ri[z]>=n||!isfinite(cv[z])){free(count);free(cursor);return -1;}if(cv[z]!=0.)count[ri[z]+1]++;}
 for(int r=0;r<n;r++){count[r+1]+=count[r];cursor[r]=count[r];}
 for(int q=0;q<n;q++)for(int64_t z=cp[q];z<cp[q+1];z++)if(cv[z]!=0.){
  int r=ri[z],at=cursor[r];if(at>count[r]&&c->columns[at-1]==q)c->values[at-1]+=cv[z];
  else{c->columns[at]=q;c->values[at]=cv[z];cursor[r]++;}
 }
 int used=0;for(int r=0;r<n;r++){
  int len=cursor[r]-count[r];c->rowptr[r]=used;
  memmove(c->columns+used,c->columns+count[r],(size_t)len*sizeof(int));
  memmove(c->values+used,c->values+count[r],(size_t)len*sizeof(double));used+=len;
 }c->rowptr[n]=used;free(count);free(cursor);c->n=n;c->source_cp=cp;return 0;
}
/* Sorted native state ownership is unchanged. Within each sector only the
 * symmetry-basis rows of locally owned states are retained, state-fast. */
static int CpuStateBind(CpuStateCoefficients*c,int ns,const int*ends,const int*inv,int wanted,int first,int last,double*data,size_t capacity){
 if(!c||!c->n||ns<1||ns>48||!ends||!inv||wanted<1||wanted>c->n||first<1||last<first-1||!data)return -1;
 free(c->positions);free(c->state_sector);free(c->state_slot);free(c->source_slot);free(c->reference);c->reference=NULL;
 int nk=last>=first?last-first+1:0;c->nk=nk;c->first=first;c->last=last;c->ns=ns;c->used=0;c->data=data;c->audit_error=0.;
 c->positions=malloc((size_t)(nk?nk:1)*sizeof(int));c->state_sector=malloc((size_t)(nk?nk:1)*sizeof(int));
 c->state_slot=malloc((size_t)(nk?nk:1)*sizeof(int));c->source_slot=malloc((size_t)c->n*sizeof(int));
 if(!c->positions||!c->state_sector||!c->state_slot||!c->source_slot)return -1;
 for(int q=0;q<c->n;q++)c->source_slot[q]=-1;
 int previous=0,states=0;
 for(int s=0;s<ns;s++){
  if(ends[s]<=previous||ends[s]>c->n)return -1;c->starts[s]=previous;c->ends[s]=ends[s];c->count[s]=0;c->position_start[s]=states;c->offset[s]=c->used;
  for(int q=previous;q<ends[s];q++)if(inv[q]>=first-1&&inv[q]<last){
   if(states>=nk)return -1;
   int local=inv[q]-(first-1),slot=c->count[s]++;c->positions[states++]=local;c->state_sector[local]=s;c->state_slot[local]=slot;c->source_slot[q]=slot;
  }
  c->used+=(size_t)(ends[s]-previous)*c->count[s];previous=ends[s];
 }c->position_start[ns]=states;
 if(previous!=c->n||states!=nk||last>wanted||c->used>capacity)return -1;return 0;
}
static int CpuStateRows(void*context,int first_row,int rows,double*out,int ld){
 CpuStateCoefficients*c=context;if(!c||!c->data||first_row<0||rows<0||first_row+rows>c->n||!out||ld<c->nk)return -1;
 for(int r=0;r<rows;r++){
  double*d=out+(size_t)r*ld;for(int k=0;k<c->nk;k++)d[k]=0.;
  int ao=first_row+r;
  for(int z=c->rowptr[ao];z<c->rowptr[ao+1];z++){
   int q=c->columns[z],s=0;while(q>=c->ends[s])s++;
   double v=c->values[z];if(v==0.||!c->count[s])continue;
   const double*source=c->data+c->offset[s]+(size_t)(q-c->starts[s])*c->count[s];
   for(int k=0;k<c->count[s];k++)d[c->positions[c->position_start[s]+k]]+=v*source[k];
  }
  for(int k=0;k<c->nk;k++){
   if(!isfinite(d[k]))return -1;
   if(c->reference){double ref=c->reference[(size_t)ao*c->nk+k];if(!isfinite(ref))return -1;double error=fabs(d[k]-ref);if(error>c->audit_error)c->audit_error=error;}
  }
 }return 0;
}
/* Explicit AO norm for the actual cut states. Never replace finite stored U
 * by an ideal orthogonality identity in the representative-tail bound. */
static int CpuStateCutNorms(CpuStateCoefficients*c,int count,const int*ids,long double*norms){
 if(!c||count<0||(count&&!ids)||(c->nk&&!norms))return -1;
 for(int k=0;k<c->nk;k++)norms[k]=NAN;
 for(int k=0;k<count;k++){
  if(ids[k]<c->first||ids[k]>c->last)continue;
  int local=ids[k]-c->first,s=c->state_sector[local],slot=c->state_slot[local];long double norm=0.;
  for(int row=0;row<c->n;row++){
   double x=0.;for(int z=c->rowptr[row];z<c->rowptr[row+1];z++){
    int q=c->columns[z];if(q<c->starts[s]||q>=c->ends[s]||c->values[z]==0.)continue;
    x+=c->values[z]*c->data[c->offset[s]+(size_t)(q-c->starts[s])*c->count[s]+slot];
   }
   if(!isfinite(x))return -1;norm+=(long double)x*(long double)x;
   if(c->reference){double error=fabs(x-c->reference[(size_t)row*c->nk+local]);if(error>c->audit_error)c->audit_error=error;}
  }if(!isfinite(norm))return -1;norms[local]=norm;
 }return 0;
}
static void CpuStateDestroy(CpuStateCoefficients*c){
 if(!c)return;free(c->rowptr);free(c->columns);free(c->values);free(c->positions);free(c->state_sector);free(c->state_slot);free(c->source_slot);free(c->reference);memset(c,0,sizeof(*c));
}
/* Route complete, unchanged retained states from reset-origin sector source
 * grids to their native state owners. Source q maps to the same global sorted
 * index as before; no occupation, degeneracy or state-count approximation. */
static int CpuStatePack(MPI_Comm world,int pr,int pc,int prow,int pcol,const CpuNativeSector*sectors,
 const double*input,const int*inv,int wanted,CpuStateCoefficients*c,int tile,size_t*temporary_bytes){
 if(!c)return -1;
 int np;MPI_Comm_size(world,&np);int valid=c&&input&&inv&&sectors&&c->data&&c->data!=input&&np==pr*pc&&tile>0&&tile<=1024;
 int*meta=calloc((size_t)(2*np+c->n),sizeof(int)),*counts=calloc((size_t)np*5+2,sizeof(int));
 if(!meta||!counts){MPI_Abort(world,96);abort();}
 int*owner_state=meta+2*np,range[2]={c->first,c->last};MPI_Allgather(range,2,MPI_INT,meta,2,MPI_INT,world);
 for(int k=0;k<c->n;k++)owner_state[k]=-1;
 for(int r=0;r<np;r++){int a=meta[2*r],b=meta[2*r+1];if(b>=a){if(a<1||b>wanted)valid=0;else for(int k=a-1;k<b;k++){if(owner_state[k]>=0)valid=0;owner_state[k]=r;}}}
 for(int k=0;k<wanted;k++)if(owner_state[k]<0)valid=0;
 int all;MPI_Allreduce(&valid,&all,1,MPI_INT,MPI_MIN,world);if(!all){free(meta);free(counts);return -1;}
 int *sc=counts,*rc=sc+np,*sd=rc+np,*rd=sd+np+1,*cursor=rd+np+1;
 CpuNativePackEntry*send=NULL,*recv=NULL;size_t sendcap=0,recvcap=0,received=0,peak=0;
 MPI_Datatype entry;MPI_Type_contiguous((int)sizeof(CpuNativePackEntry),MPI_BYTE,&entry);MPI_Type_commit(&entry);int bad=0;
 for(int begin=0;begin<c->n;begin+=tile){
  int end=c->n-begin<tile?c->n:begin+tile;memset(sc,0,(size_t)np*sizeof(int));
  for(int s=0;s<c->ns;s++){
   const CpuNativeSector*t=sectors+s;for(int j=0;j<t->nc;j++){
    int q=t->ofs+CpuNativePackGlobal(j,t->nb,pcol,pc);if(q>=t->ofs+t->n||q<begin||q>=end||inv[q]>=wanted)continue;
    int dest=owner_state[inv[q]];if((int64_t)sc[dest]+t->nr>INT_MAX){MPI_Abort(world,96);abort();}sc[dest]+=t->nr;
   }
  }
  MPI_Alltoall(sc,1,MPI_INT,rc,1,MPI_INT,world);sd[0]=rd[0]=0;
  for(int r=0;r<np;r++){if((int64_t)sd[r]+sc[r]>INT_MAX||(int64_t)rd[r]+rc[r]>INT_MAX){MPI_Abort(world,96);abort();}sd[r+1]=sd[r]+sc[r];rd[r+1]=rd[r]+rc[r];cursor[r]=sd[r];}
  if((size_t)sd[np]>sendcap){free(send);sendcap=sd[np];send=malloc((sendcap?sendcap:1)*sizeof(*send));}
  if((size_t)rd[np]>recvcap){free(recv);recvcap=rd[np];recv=malloc((recvcap?recvcap:1)*sizeof(*recv));}
  if(!send)send=malloc(sizeof(*send));if(!recv)recv=malloc(sizeof(*recv));if(!send||!recv){MPI_Abort(world,96);abort();}
  size_t bytes=(sendcap+recvcap)*sizeof(*send)+(size_t)(np*7+2+c->n)*sizeof(int);if(bytes>peak)peak=bytes;
  for(int s=0;s<c->ns;s++){
   const CpuNativeSector*t=sectors+s;for(int j=0;j<t->nc;j++){
    int q=t->ofs+CpuNativePackGlobal(j,t->nb,pcol,pc);if(q>=t->ofs+t->n||q<begin||q>=end||inv[q]>=wanted)continue;
    int dest=owner_state[inv[q]];
    for(int i=0;i<t->nr;i++){int row=CpuNativePackGlobal(i,t->nb,prow,pr);if(row>=t->n)continue;send[cursor[dest]++]=(CpuNativePackEntry){row,q,input[t->offset+(size_t)j*t->nr+i]};}
   }
  }
  for(int r=0;r<np;r++)if(cursor[r]!=sd[r+1])bad=1;
  MPI_Alltoallv(send,sc,sd,entry,recv,rc,rd,entry,world);received+=(size_t)rd[np];
  for(int k=0;k<rd[np];k++){
   CpuNativePackEntry*v=recv+k;int s=0;while(s<c->ns&&v->col>=c->ends[s])s++;
   if(s==c->ns||v->col<c->starts[s]||v->row<0||v->row>=c->ends[s]-c->starts[s]||c->source_slot[v->col]<0||!isfinite(v->value)){bad=1;continue;}
   c->data[c->offset[s]+(size_t)v->row*c->count[s]+c->source_slot[v->col]]=v->value;
  }
 }
 if(received!=c->used)bad=1;MPI_Type_free(&entry);free(send);free(recv);free(meta);free(counts);
 MPI_Allreduce(&bad,&all,1,MPI_INT,MPI_MAX,world);if(temporary_bytes)*temporary_bytes=peak;return all?-1:0;
}
#endif
