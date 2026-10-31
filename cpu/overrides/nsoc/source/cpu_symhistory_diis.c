#include "cpu_symhistory_diis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

struct CpuSymHistoryDIIS {
  size_t length;
  int history, slots, last_step;
  MPI_Comm comm;
  double *weights, **density, **residual;
  double *gram_local, *gram_global, *a_data, *ia_data, **a, **ia, *alpha;
};

void CpuSymHistoryDIIS_Free(CpuSymHistoryDIIS *p) {
  int s;
  if(!p)return;
  if(p->density)for(s=0;s<p->slots;s++)free(p->density[s]);
  if(p->residual)for(s=0;s<p->slots;s++)free(p->residual[s]);
  free(p->density);free(p->residual);free(p->weights);
  free(p->gram_local);free(p->gram_global);free(p->a_data);free(p->ia_data);
  free(p->a);free(p->ia);free(p->alpha);free(p);
}

CpuSymHistoryDIIS *CpuSymHistoryDIIS_Create(size_t n,const double *w,int h,MPI_Comm comm) {
  int s;size_t k,pairs,mat;CpuSymHistoryDIIS *p;
  if(h<2||h>1024||(n&&!w)||n>((size_t)-1)/sizeof(double))return NULL;
  for(k=0;k<n;k++)if(!isfinite(w[k])||w[k]<=0)return NULL;
  p=calloc(1,sizeof(*p));if(!p)return NULL;
  p->length=n;p->history=h;p->slots=h+2;p->comm=comm;
  p->density=calloc(p->slots,sizeof(double*));p->residual=calloc(p->slots,sizeof(double*));
  p->weights=malloc((n?n:1)*sizeof(double));pairs=(size_t)h*(h+1)/2;mat=(size_t)p->slots*p->slots;
  p->gram_local=calloc(pairs,sizeof(double));p->gram_global=calloc(pairs,sizeof(double));
  p->a_data=calloc(mat,sizeof(double));p->ia_data=calloc(mat,sizeof(double));
  p->a=calloc(p->slots,sizeof(double*));p->ia=calloc(p->slots,sizeof(double*));p->alpha=calloc(p->slots,sizeof(double));
  if(!p->density||!p->residual||!p->weights||!p->gram_local||!p->gram_global||!p->a_data||!p->ia_data||!p->a||!p->ia||!p->alpha)goto bad;
  if(n)memcpy(p->weights,w,n*sizeof(double));
  for(s=0;s<p->slots;s++){
    p->density[s]=calloc(n?n:1,sizeof(double));p->residual[s]=calloc(n?n:1,sizeof(double));
    p->a[s]=p->a_data+(size_t)s*p->slots;p->ia[s]=p->ia_data+(size_t)s*p->slots;
    if(!p->density[s]||!p->residual[s])goto bad;
  }
  return p;
bad:CpuSymHistoryDIIS_Free(p);return NULL;
}
size_t CpuSymHistoryDIIS_Length(const CpuSymHistoryDIIS *p){return p?p->length:0;}
double *CpuSymHistoryDIIS_Density(CpuSymHistoryDIIS *p,int s){return p&&s>=0&&s<p->slots?p->density[s]:NULL;}
double *CpuSymHistoryDIIS_Residual(CpuSymHistoryDIIS *p,int s){return p&&s>=0&&s<p->slots?p->residual[s]:NULL;}
int CpuSymHistoryDIIS_Set(CpuSymHistoryDIIS *p,int r,int s,const double *v) {
  size_t k;double *dest;
  if(!p||(r!=0&&r!=1)||s<0||s>=p->slots||(p->length&&!v))return-1;
  for(k=0;k<p->length;k++)if(!isfinite(v[k]))return-3;
  dest=r?p->residual[s]:p->density[s];
  if(p->length)memcpy(dest,v,p->length*sizeof(double));return 0;
}

int CpuSymHistoryDIIS_Step(CpuSymHistoryDIIS *p,int step,double mix,CpuSymHistoryInverse inverse,void *ctx,CpuSymHistoryDIISInfo *info) {
  int m,a,b,pair,nmix,nslide,npair,bad=0,anybad=0;size_t k,bytes;double coef,localnorm=0,globalnorm=0;
  CpuSymHistoryDIISInfo out;
  if(!p||!info||step<2||(!p->last_step&&step>3)||(p->last_step&&step!=p->last_step+1)||!isfinite(mix)||!inverse)return-1;
  memset(&out,0,sizeof(out));out.effective_step=step;bytes=p->length*sizeof(double);
  /* Set() checks imported vectors; direct pointer writers are supported and
   * checked collectively so a bad owner cannot leave peers inside Allreduce. */
  for(k=0;k<p->length;k++)if(!isfinite(p->density[0][k])||!isfinite(p->density[1][k]))bad=1;
  MPI_Allreduce(&bad,&anybad,1,MPI_INT,MPI_MAX,p->comm);if(anybad)return-3;
  if(step==2){
    for(k=0;k<p->length;k++){
      double r=p->density[0][k]-p->density[1][k];p->residual[2][k]=r;
      localnorm+=p->weights[k]*r*r;
      p->density[2][k]=p->density[1][k];
      p->density[0][k]=mix*p->density[0][k]+(1.-mix)*p->density[1][k];
      p->density[1][k]=p->density[0][k];
    }
    MPI_Allreduce(&localnorm,&globalnorm,1,MPI_DOUBLE,MPI_SUM,p->comm);
    out.norm_current=globalnorm;out.mixing_coefficient=mix;
    p->last_step=step;*info=out;return 0;
  }
  for(k=0;k<p->length;k++)p->residual[1][k]=p->density[0][k]-p->density[1][k];
  if(step-1<p->history){nmix=step-1;nslide=nmix+1;}else{nmix=p->history;nslide=nmix;}
  npair=nmix*(nmix+1)/2;memset(p->gram_local,0,(size_t)npair*sizeof(double));
  /* Same full Frobenius metric in orthogonal block-orbit coordinates:
   * each stored AO entry has its directed block-orbit multiplicity weight. */
  for(k=0;k<p->length;k++){
    pair=0;
    for(a=1;a<=nmix;a++)for(b=a;b<=nmix;b++)p->gram_local[pair++]+=p->weights[k]*p->residual[a][k]*p->residual[b][k];
  }
  MPI_Allreduce(p->gram_local,p->gram_global,npair,MPI_DOUBLE,MPI_SUM,p->comm);
  pair=0;for(a=0;a<nmix;a++)for(b=a;b<nmix;b++)p->a[a][b]=p->a[b][a]=p->gram_global[pair++];
  out.norm_current=p->a[0][0];
  for(a=0;a<nmix;a++)p->a[a][nmix]=p->a[nmix][a]=-1.;p->a[nmix][nmix]=0.;
  inverse(nmix,p->a,p->ia,ctx);
  for(m=1;m<=nmix;m++){p->alpha[m]=-p->ia[m-1][nmix];if(!isfinite(p->alpha[m]))out.fallback=1;}
  if(out.fallback){for(m=1;m<=nmix;m++)p->alpha[m]=0.;p->alpha[1]=.1;p->alpha[2]=.9;}
  if(out.norm_current>=1e-1)coef=.5;
  else if(out.norm_current>=1e-2)coef=.6;
  else if(out.norm_current>=1e-3)coef=.7;
  else if(out.norm_current>=1e-4)coef=.8;
  else coef=1.;
  for(k=0;k<p->length;k++){
    double sum_r=0,sum_d=0;
    for(m=1;m<=nmix;m++){sum_r+=p->alpha[m]*p->residual[m][k];sum_d+=p->alpha[m]*p->density[m][k];}
    localnorm+=p->weights[k]*sum_r*sum_r;p->density[0][k]=sum_d+coef*sum_r;
  }
  MPI_Allreduce(&localnorm,&globalnorm,1,MPI_DOUBLE,MPI_SUM,p->comm);out.norm_optimized=globalnorm;
  /* Keep slot0 stable. Native density shift includes mixed D0->D1. The
   * residual's live R1 stays intact and is duplicated into R2, exactly as
   * the original copy loops (rather than recycling R1 before its next read). */
  {double *recycled=p->density[nslide];for(m=nslide;m>1;m--)p->density[m]=p->density[m-1];p->density[1]=recycled;if(bytes)memcpy(p->density[1],p->density[0],bytes);}
  {double *recycled=p->residual[nslide];for(m=nslide;m>2;m--)p->residual[m]=p->residual[m-1];p->residual[2]=recycled;if(bytes)memcpy(p->residual[2],p->residual[1],bytes);}
  out.num_mix=nmix;out.num_slide=nslide;out.mixing_coefficient=coef;p->last_step=step;*info=out;return 0;
}
