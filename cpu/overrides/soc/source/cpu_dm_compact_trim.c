#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <mkl.h>
#include "cpu_dm_compact_trim.h"
#include <stdio.h>
#include <time.h>
extern int CpuC3_CompactPack(int row0,int rows,int band0,int bands,void*out,int ld);

static int dm_trim_enabled=0,dm_profile_enabled=0;
static CpuDMCompactProfile dm_profile_last;
static double dm_trim_now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1.e-9*t.tv_nsec;}
void CpuDMCompactSetTrim(int enabled){if(enabled!=0&&enabled!=1){fprintf(stderr,"CpuDMCompactSetTrim requires0 or1\n");abort();}dm_trim_enabled=enabled;}
void CpuDMCompactSetProfile(int enabled){if(enabled!=0&&enabled!=1){fprintf(stderr,"CpuDMCompactSetProfile requires0 or1\n");abort();}dm_profile_enabled=enabled;}
void CpuDMCompactGetProfile(CpuDMCompactProfile *out){if(out)*out=dm_profile_last;}

int cpu_dm_blocks_compact(int n,int bands,const void *vectors,const double*w,
 int count,const int*blocks,int entries,int edm,int tile,double*out)
{
 int p,c,first,last,m,nc,maxm=0,maxc=0,i,j,k,t,nb,pair,base,b,offset;
 int components=edm?4:6,rc=0;
 double kernel_start=0.,tick=0.;
 if(dm_profile_enabled){memset(&dm_profile_last,0,sizeof(dm_profile_last));dm_profile_last.logical_bands=bands;kernel_start=dm_trim_now();}
 MKL_Complex16 *A=NULL,*B=NULL,*G=NULL;
 const MKL_Complex16 one={1.,0.};
 (void)vectors;
 if(n<0||n>INT_MAX/2||bands<0||count<0||entries<0||tile<=0||!(edm==0||edm==1))return -1;
 if((count&&!blocks)||(entries&&!out)||(bands&&!w))return -1;
 for(c=0;c<count;c++){
  const int*x=blocks+5*c;
  if(x[0]<0||x[1]<=0||x[2]<0||x[3]<=0||x[4]<0||
     (int64_t)x[0]+x[1]>n||(int64_t)x[2]+x[3]>n||
     (int64_t)x[4]+(int64_t)x[1]*x[3]>entries)return -1;
 }
 /* Exact tail zeros only: preserve the existing occupation filter, signed
  * energy weights, tiny nonzero values and internal zero panels unchanged. */
 if(dm_trim_enabled)while(bands>0&&w[bands-1]==0.)bands--;
 if(dm_profile_enabled){dm_profile_last.effective_bands=bands;for(k=0;k<bands;k++)if(w[k]!=0.)dm_profile_last.nonzero_weights++;}
 if(entries)memset(out,0,(size_t)components*entries*sizeof(double));
 if(!bands||!count){if(dm_profile_enabled)dm_profile_last.kernel_seconds=dm_trim_now()-kernel_start;return 0;}
 for(k=0;k<bands&&w[k]==0.;k++);
 if(k==bands){if(dm_profile_enabled)dm_profile_last.kernel_seconds=dm_trim_now()-kernel_start;return 0;}
 for(first=0;first<count;first=last){
  m=blocks[5*first+1];nc=0;last=first;
  while(last<count&&blocks[5*last]==blocks[5*first]&&blocks[5*last+1]==m){
   if((int64_t)nc+blocks[5*last+3]>2147483647)return -1;
   nc+=blocks[5*last+3];last++;
  }
  if(m>maxm)maxm=m;if(nc>maxc)maxc=nc;
 }
 if(tile>bands)tile=bands;
 if((size_t)maxm>SIZE_MAX/sizeof(*A)/(size_t)tile||
    (size_t)maxc>SIZE_MAX/sizeof(*B)/(size_t)tile||
    (size_t)maxm>SIZE_MAX/sizeof(*G)/(size_t)maxc)return -1;
 A=malloc((size_t)maxm*tile*sizeof(*A));B=malloc((size_t)maxc*tile*sizeof(*B));G=malloc((size_t)maxm*maxc*sizeof(*G));
 if(!A||!B||!G){rc=-2;goto done;}
 if(dm_profile_enabled)dm_profile_last.scratch_bytes=sizeof(*A)*((size_t)maxm*tile+(size_t)maxc*tile+(size_t)maxm*maxc);
 for(first=0;first<count;first=last){
  m=blocks[5*first+1];nc=0;last=first;
  while(last<count&&blocks[5*last]==blocks[5*first]&&blocks[5*last+1]==m){nc+=blocks[5*last+3];last++;}
  for(pair=0;pair<3;pair++){
   int sa=pair==1?n:0,sb=pair==0?0:n;
   memset(G,0,(size_t)m*nc*sizeof(*G));
   for(k=0;k<bands;k+=tile){
    nb=bands-k<tile?bands-k:tile;
    /* Compact storage includes all logical states. Do not unpack completely
     * empty occupation tiles, including states beyond the requested MaxN. */
    for(t=0;t<nb&&w[k+t]==0.;t++);
    if(t==nb)continue;
    if(dm_profile_enabled){tick=dm_trim_now();dm_profile_last.pack_calls++;}
    if(CpuC3_CompactPack(sa+blocks[5*first],m,k,nb,A,m)){rc=-3;goto done;}
    base=0;
    for(b=first;b<last;b++){
     if(dm_profile_enabled)dm_profile_last.pack_calls++;
     if(CpuC3_CompactPack(sb+blocks[5*b+2],blocks[5*b+3],k,nb,B+base,nc)){rc=-3;goto done;}
     base+=blocks[5*b+3];
    }
    for(t=0;t<nb;t++){
     if(w[k+t]!=0.){for(i=0;i<m;i++){A[(size_t)t*m+i].real*=w[k+t];A[(size_t)t*m+i].imag*=w[k+t];}}
     else{memset(A+(size_t)t*m,0,(size_t)m*sizeof(*A));memset(B+(size_t)t*nc,0,(size_t)nc*sizeof(*B));}
    }
    if(dm_profile_enabled){dm_profile_last.pack_seconds+=dm_trim_now()-tick;tick=dm_trim_now();dm_profile_last.gemm_calls++;}
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,m,nc,nb,&one,A,m,B,nc,&one,G,m);
    if(dm_profile_enabled)dm_profile_last.gemm_seconds+=dm_trim_now()-tick;
   }
   if(dm_profile_enabled)tick=dm_trim_now();
   base=0;
   for(b=first;b<last;b++){
    int q=blocks[5*b+3];offset=blocks[5*b+4];
    for(i=0;i<m;i++)for(j=0;j<q;j++){
     MKL_Complex16 z=G[i+(size_t)(base+j)*m];p=offset+i*q+j;
     if(pair==0){out[p]=z.real;if(!edm)out[(size_t)4*entries+p]=-z.imag;}
     if(pair==1){out[(size_t)entries+p]=z.real;if(!edm)out[(size_t)5*entries+p]=-z.imag;}
     if(pair==2){out[(size_t)2*entries+p]=z.real;out[(size_t)3*entries+p]=-z.imag;}
    }
    base+=q;
   }
   if(dm_profile_enabled)dm_profile_last.scatter_seconds+=dm_trim_now()-tick;
  }
 }
done:
 free(A);free(B);free(G);
 if(dm_profile_enabled)dm_profile_last.kernel_seconds=dm_trim_now()-kernel_start;
 return rc;
}
