/* Generic GSOC independent-sector coefficient store.
 * Column ownership j % np keeps each independent vector once across MPI.
 * Physical/TR logical states refer to that vector; no dense physical EVec is
 * built. Pack is local and preserves every logical band's individual weight. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include "mpi.h"
#include "openmx_common.h"
#include "cpu_gsoc_compact.h"

typedef struct {int n,ofs,dtype,partner,child,nowned;void *q;} CSector;
typedef struct {int sorted,sector,column;} CBand;
typedef struct {double e;int source;} CEigen;
static struct {
 int built,ns,nsec,np,rank,maxn,nlocal;
 MPI_Comm comm;
 CSector *sec;
 CBand *band;
 int *rowptr,*column,*global_owner,*global_slot;
 dcomplex *value;
 size_t bytes,qbytes;
} C;
static void cdie(const char *s){int rank=0;MPI_Comm_rank(C.built?C.comm:MPI_COMM_WORLD,&rank);fprintf(stderr,"CPU_GSOC_COMPACT_ERROR rank=%d %s\n",rank,s);fflush(stderr);MPI_Abort(C.built?C.comm:MPI_COMM_WORLD,95);abort();}
static void *cmem(size_t n,size_t bytes){if(bytes&&n>SIZE_MAX/bytes)cdie("allocation capacity overflow");void *p=calloc(n?n:1,bytes);if(!p)cdie("allocation failed");return p;}
static int cglobal(int i,int nb,int p,int np){return(i/nb*np+p)*nb+i%nb;}
static dcomplex cz(double r,double i){dcomplex z;z.r=r;z.i=i;return z;}
static dcomplex cmul(dcomplex a,dcomplex b){return cz(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r);}
static int cecomp(const void *aa,const void *bb){const CEigen *a=aa,*b=bb;if(a->e<b->e)return-1;if(a->e>b->e)return 1;return(a->source>b->source)-(a->source<b->source);}
int CpuGSOC_CompactEnabled(void){const char *s=getenv("OPENMX_CPU_GSOC_COMPACT");if(!s||!strcmp(s,"0"))return 0;if(strcmp(s,"1"))cdie("COMPACT must be exactly 0/1");return 1;}
int CpuGSOC_CompactLocalBands(void){return C.built?C.nlocal:0;}
int CpuGSOC_CompactSortedBand(int b){if(!C.built||b<0||b>=C.nlocal)cdie("logical local band out of range");return C.band[b].sorted+1;}
size_t CpuGSOC_CompactBytes(void){return C.bytes;}

/* Deterministic wire ordering avoids element coordinates. Each sender visits
 * local columns then local rows. The receiver derives those same indices from
 * the BLACS rank; only real/complex FP64 values cross the network. */
static void ctransfer(CpuGSOCCompactSector *in,CSector *out){
 int p=C.np,fields=in->dtype==1?1:2,n=in->n;
 out->nowned=C.rank<n?(n-1-C.rank)/p+1:0;
 size_t item=in->dtype==1?sizeof(double):sizeof(dcomplex);
 size_t count=(size_t)out->nowned*n;
 out->q=cmem(count,item);
 if(count>SIZE_MAX/item||C.bytes>SIZE_MAX-count*item)cdie("coefficient byte count overflow");C.bytes+=count*item;C.qbytes+=count*item;
 int *ct=cmem(p,sizeof(int)),*rc=cmem(p,sizeof(int)),*dp=cmem(p,sizeof(int)),*rd=cmem(p,sizeof(int)),*pos=cmem(p,sizeof(int));
 for(int b0=0;b0<n;b0+=16){
  int b1=n-b0<16?n:b0+16;memset(ct,0,(size_t)p*sizeof(int));
  for(int j=0;j<in->nc;j++){int b=cglobal(j,in->nb,in->pcol,in->npcols);if(b<b0||b>=b1)continue;int dest=b%p;
   for(int i=0;i<in->nr;i++)if(cglobal(i,in->nb,in->prow,in->nprows)<n){if(ct[dest]>INT_MAX-fields)cdie("column transfer MPI count overflow");ct[dest]+=fields;}}
  MPI_Alltoall(ct,1,MPI_INT,rc,1,MPI_INT,C.comm);dp[0]=rd[0]=0;
  for(int r=1;r<p;r++){if(dp[r-1]>INT_MAX-ct[r-1]||rd[r-1]>INT_MAX-rc[r-1])cdie("column transfer MPI displacement overflow");dp[r]=dp[r-1]+ct[r-1];rd[r]=rd[r-1]+rc[r-1];}
  if(dp[p-1]>INT_MAX-ct[p-1]||rd[p-1]>INT_MAX-rc[p-1])cdie("column transfer MPI total overflow");
  int sendn=dp[p-1]+ct[p-1],recvn=rd[p-1]+rc[p-1];double *send=cmem(sendn,sizeof(double)),*recv=cmem(recvn,sizeof(double));memcpy(pos,dp,(size_t)p*sizeof(int));
  for(int j=0;j<in->nc;j++){int b=cglobal(j,in->nb,in->pcol,in->npcols);if(b<b0||b>=b1)continue;int dest=b%p;
   for(int i=0;i<in->nr;i++)if(cglobal(i,in->nb,in->prow,in->nprows)<n){size_t at=(size_t)j*in->nr+i;int to=pos[dest];if(to>dp[dest]+ct[dest]-fields)cdie("column pack capacity");
    if(fields==1)send[to++]=((double*)in->Q)[at];else{dcomplex z=((dcomplex*)in->Q)[at];send[to++]=z.r;send[to++]=z.i;}pos[dest]=to;}}
  MPI_Alltoallv(send,ct,dp,MPI_DOUBLE,recv,rc,rd,MPI_DOUBLE,C.comm);
  for(int source=0;source<p;source++){
   int pr=source/in->npcols,pc=source%in->npcols,at=rd[source];
   for(int b=b0;b<b1;b++)if(b%p==C.rank&&(b/in->nb)%in->npcols==pc){
    for(int row=0;row<n;row++)if((row/in->nb)%in->nprows==pr){size_t to=(size_t)(b/p)*n+row;if(to>=count||at>rd[source]+rc[source]-fields)cdie("column unpack capacity");
     if(fields==1)((double*)out->q)[to]=recv[at++];else{((dcomplex*)out->q)[to]=cz(recv[at],recv[at+1]);at+=2;}}}
   if(at!=rd[source]+rc[source])cdie("column transfer sender/receiver ordering mismatch");}
  free(send);free(recv);
 }
 free(ct);free(rc);free(dp);free(rd);free(pos);free(in->Q);in->Q=NULL;
}

void CpuGSOC_CompactBuild(MPI_Comm comm,int ns,int nsec,CpuGSOCCompactSector *sec,const int64_t *cp,const int *ri,const dcomplex *cv,const double *eval,int maxn,double *ko){
 if(C.built)cdie("provider already built");if(!CpuGSOC_CompactEnabled())cdie("disabled provider build");
 if(ns<1||nsec<1||nsec>48||maxn<1||maxn>ns||!sec||!cp||!ri||!cv||!eval||!ko)cdie("invalid build arguments");
 C.comm=comm;C.built=1;C.ns=ns;C.nsec=nsec;C.maxn=maxn;MPI_Comm_rank(comm,&C.rank);MPI_Comm_size(comm,&C.np);
 C.sec=cmem(nsec,sizeof(CSector));C.rowptr=cmem((size_t)ns+1,sizeof(int));
 if(cp[0]!=0||cp[ns]<ns||cp[ns]>INT_MAX)cdie("invalid sparse basis extent");
 int nnz=(int)cp[ns];C.column=cmem(nnz,sizeof(int));C.value=cmem(nnz,sizeof(dcomplex));
 C.bytes=(size_t)nsec*sizeof(CSector)+((size_t)ns+1)*sizeof(int)+(size_t)nnz*(sizeof(int)+sizeof(dcomplex));
 for(int j=0;j<ns;j++){if(cp[j]<0||cp[j]>=cp[j+1]||cp[j+1]>nnz)cdie("invalid sparse basis pointer");for(int64_t a=cp[j];a<cp[j+1];a++){int row=ri[a];if(row<0||row>=ns||!isfinite(cv[a].r)||!isfinite(cv[a].i))cdie("invalid sparse basis element");C.rowptr[row+1]++;}}
 for(int i=0;i<ns;i++)C.rowptr[i+1]+=C.rowptr[i];int *next=cmem(ns,sizeof(int));memcpy(next,C.rowptr,(size_t)ns*sizeof(int));
 for(int j=0;j<ns;j++)for(int64_t a=cp[j];a<cp[j+1];a++){int at=next[ri[a]]++;C.column[at]=j;C.value[at]=cv[a];}free(next);
 CEigen *sort=cmem(ns,sizeof(CEigen));int *inv=cmem(ns,sizeof(int)),*which=cmem(ns,sizeof(int));
 int end=0;
 for(int k=0;k<nsec;k++){
  CpuGSOCCompactSector *s=sec+k;CSector *t=C.sec+k;
  if(s->ofs!=end||s->n<1||s->n>ns-end||(s->dtype!=1&&s->dtype!=2)||s->partner<-1||s->partner>=k||s->nprows*s->npcols!=C.np||s->prow*s->npcols+s->pcol!=C.rank||s->nb<1)cdie("invalid sector geometry");
  t->n=s->n;t->ofs=s->ofs;t->dtype=s->dtype;t->partner=s->partner;t->child=s->child;
  if(s->partner>=0){CSector *p=C.sec+s->partner;if(p->partner>=0||p->n!=s->n||p->dtype!=s->dtype||p->child!=k)cdie("invalid TR sector graph");}
  for(int j=0;j<s->n;j++)which[end+j]=k;end+=s->n;
 }
 if(end!=ns)cdie("sector sizes do not span spinor space");
 for(int j=0;j<ns;j++){if(!isfinite(eval[j]))cdie("nonfinite eigenvalue");sort[j].e=eval[j];sort[j].source=j;}
 qsort(sort,ns,sizeof(CEigen),cecomp);for(int j=0;j<ns;j++){ko[j+1]=sort[j].e;inv[sort[j].source]=j;}
 C.global_owner=cmem(maxn,sizeof(int));C.global_slot=cmem(maxn,sizeof(int));int *owner_count=cmem(C.np,sizeof(int));
 for(int b=0;b<maxn;b++){int source=sort[b].source,k=which[source],col=source-C.sec[k].ofs,owner=col%C.np;C.global_owner[b]=owner;C.global_slot[b]=owner_count[owner]++;}
 C.nlocal=owner_count[C.rank];C.band=cmem(C.nlocal,sizeof(CBand));C.bytes+=(size_t)maxn*2*sizeof(int)+(size_t)C.nlocal*sizeof(CBand);
 int local=0;for(int b=0;b<maxn;b++)if(C.global_owner[b]==C.rank){int source=sort[b].source,k=which[source];C.band[local].sorted=b;C.band[local].sector=k;C.band[local].column=source-C.sec[k].ofs;local++;}
 for(int k=0;k<nsec;k++)if(sec[k].partner<0){if(!sec[k].Q)cdie("independent Q missing");ctransfer(sec+k,C.sec+k);}
 unsigned long long own[3]={C.bytes,C.qbytes,C.bytes-C.qbytes},total[3]={0};MPI_Reduce(own,total,3,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,comm);
 if(!C.rank){printf("CPU_GSOC_COMPACT_STORE ns=%d requested=%d local_state_layout=cyclic_independent_column independent_Q_bytes_plus_metadata=%llu Q_bytes_global=%llu metadata_bytes_global=%llu physical_EVec_bytes=0 partner_Q_bytes=0\n",ns,maxn,total[0],total[1],total[2]);fflush(stdout);}
 free(sort);free(inv);free(which);free(owner_count);
}

int CpuGSOC_CompactPack(int row0,int rows,int band0,int bands,void *output,int ld){
 if(!C.built||row0<0||rows<0||row0>C.ns-rows||band0<0||bands<0||band0>C.nlocal-bands||ld<rows||(!output&&rows&&bands))return-1;
 dcomplex *out=output;
 for(int b=0;b<bands;b++){
  CBand *v=C.band+band0+b;CSector *logical=C.sec+v->sector,*source=logical->partner<0?logical:C.sec+logical->partner;
  if(v->column%C.np!=C.rank)return-1;size_t qbase=(size_t)(v->column/C.np)*source->n;
  for(int i=0;i<rows;i++){
   int row=row0+i;dcomplex sum=cz(0,0);
   for(int a=C.rowptr[row];a<C.rowptr[row+1];a++){
    int qrow=C.column[a]-logical->ofs;if(qrow<0||qrow>=logical->n)continue;
    dcomplex q=source->dtype==1?cz(((double*)source->q)[qbase+qrow],0):((dcomplex*)source->q)[qbase+qrow];if(logical->partner>=0)q.i=-q.i;
    dcomplex value=cmul(C.value[a],q);sum.r+=value.r;sum.i+=value.i;
   }
   out[(size_t)b*ld+i]=sum;
  }
 }
 return 0;
}
int CpuGSOC_CompactGetPhysical(int band1,dcomplex *out){
 if(!C.built||band1<1||band1>C.maxn||!out)return-1;int b=band1-1,owner=C.global_owner[b],bad=0;
 if(C.rank==owner)bad=CpuGSOC_CompactPack(0,C.ns,C.global_slot[b],1,out,C.ns);
 MPI_Bcast(&bad,1,MPI_INT,owner,C.comm);if(bad)return-1;MPI_Bcast(out,2*C.ns,MPI_DOUBLE,owner,C.comm);return 0;
}
void CpuGSOC_CompactRelease(void){
 for(int k=0;k<C.nsec;k++)free(C.sec[k].q);free(C.sec);free(C.band);free(C.rowptr);free(C.column);free(C.value);free(C.global_owner);free(C.global_slot);memset(&C,0,sizeof(C));
}
